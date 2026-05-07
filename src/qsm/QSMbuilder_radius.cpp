/**
 * @file QSMbuilder_radius.cpp
 * Project: Arbor
 * 
 * Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <map>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <iomanip>

#include "QSMbuilder.h"
#include "PointCloud.h"
#include "nanoflann.h"
#include "fitting.h"

namespace arbor::qsm {

class SimpleAdaptor
{
public:
  struct Point { double x, y, z; int id; };
  std::vector<Point> points;
  inline size_t kdtree_get_point_count() const { return points.size(); }
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const
  {
    if (dim == 0) return points[idx].x;
    if (dim == 1) return points[idx].y;
    return points[idx].z;
  }
  template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
  inline size_t point_count() const { return points.size(); }
  inline size_t size() const { return points.size(); }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, SimpleAdaptor>, SimpleAdaptor, 3> CentroidKDTree;

// Griese, N., Ritzert, M. & Nölke, N. A large dataset of labelled single tree point clouds, QSMs and
// tree graphs. Sci Data 12, 1953 (2025). https://doi.org/10.1038/s41597-025-06421-7
class Allometry
{
public:
  static double H_vs_DBH(double dbh)
  {
    return 36.03 * std::pow(1.0 - std::exp(-0.05 * dbh), 1.1);
  }

  static double DBH_vs_H(double H)
  {
    double DBH_cm;
    if (H < 25.0)
    {
      double ratio = std::clamp(H / 36.03, 0.0, 1.0 - 1e-12);
      DBH_cm = -1.0 / 0.05 * std::log(1.0 - std::pow(ratio, 1.0 / 1.1));
    }
    else {
      DBH_cm = 4.0 * H - 75.0;
    }
    return DBH_cm / 100.0;
  }
};

void QSMbuilder::construct_radii(const PointCloud& tree, double tip_radius)
{
  double H = 0.0;
  for (size_t i = 0; i < tree.size(); ++i)
  {
    double z = tree.get_hag(i);
    if (z > H) H = z;
  }

  // Estimation of an expected radius based on broad allometry
  double R0 = Allometry::DBH_vs_H(H) / 2.0;

  // TODO: control parameter
  // If the estimated radius is too small we don't even try to measure the tree
  if (R0 < 0.03)
  {
    conic_allometry(R0, tip_radius);

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2);
    oss << "[WARN 1] H = " << H
        << std::fixed << std::setprecision(1)
        << " m: estimated DBH = " << (2 * R0 * 100)
        << " cm. Too small to be measured. "
        << "The QSM is the result of pure allometry without actual measurements.";

    std::string msg = oss.str();

    graph.messages.push_back(msg);
    ServiceLocator::logger()("\033[33m" + msg + "\033[0m");
    return;
  }

  ServiceLocator::logger()("Pre-allometry");
  conic_allometry(2.0 * R0, tip_radius);

  ServiceLocator::logger()("Measuring diameters");
  measure_radii(tree, 180.0, 0.2, 0.3, 0.03);

  compute_architecture(true);

  ServiceLocator::logger()("Polynomial fitting");
  polynomial_fitting(tip_radius);


  // Check if main axis has valid measurements
  bool has_na = false;
  for (auto& [eid, einfo] : graph.edges())
  {
    if (einfo.data.axis_ID == 1 && einfo.data.radius == RADIUS_UNSET)
    {
      has_na = true;
      break;
    }
  }

  if (has_na)
  {
    std::string msg = "[WARN 2] Not a single valid measure for this tree. The QSM is a pure reconstruction based on allometry";
    ServiceLocator::logger()("\033[33m" + msg + "\033[0m");
    graph.messages.push_back(msg);
    conic_allometry(R0, tip_radius);
    return;
  }


  // Check root radius
  double Rroot;
  for (auto& [eid, einfo] : graph.edges())
  {
    if (einfo.data.axis_ID == 1 && einfo.data.parent_ID == 0)
    {
      Rroot = einfo.data.radius;
      break;
    }
  }

  /*if (Rroot > 3*R0)
  {
    std::ostringstream oss;
    oss << "[WARN 3] Measured root diameter is "
        << static_cast<int>(Rroot / R0)
        << " times greater than the expected DBH ("
        << std::fixed << std::setprecision(1) << 2 * Rroot
        << " vs. "
        << std::fixed << std::setprecision(2) << 2 * R0
        << ").  Fall back to pure reconstruction based on allometry";

    std::string msg = oss.str();
    ServiceLocator::logger()("\033[33m" + msg + "\033[0m");
    graph.messages.push_back(msg);

    conic_allometry(2.0 * R0, tip_radius);
    return;
  }*/

  ServiceLocator::logger()("Reconstruction");
  reconstruct_missing_radii(tip_radius);
}

void QSMbuilder::measure_radii(const PointCloud& tree, float sarc, float sins, float sinl, float srmeas)
{
  if (graph.edge_count() == 0) throw std::runtime_error("Internal error: no cylinder in this QSM. Please report.");

  // Prepare centroids for KD-Tree
  SimpleAdaptor centroids_cloud;
  centroids_cloud.points.reserve(graph.edge_count());

  std::vector<int> index_to_eid;
  index_to_eid.reserve(graph.edge_count());

  for (auto& [eid, einfo] : graph.edges())
  {
    einfo.data.radius = RADIUS_UNSET;

    const QSMNode& src = graph.node(einfo.source);
    const QSMNode& tgt = graph.node(einfo.target);

    double cx = (src.x + tgt.x) / 2.0;
    double cy = (src.y + tgt.y) / 2.0;
    double cz = (src.z + tgt.z) / 2.0;

    centroids_cloud.points.push_back({cx, cy, cz, eid});
    index_to_eid.push_back(eid);
  }

  CentroidKDTree index(3, centroids_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  // Assign tree points to nearest edge (cylinder)
  std::vector<std::vector<size_t>> points_per_edge(centroids_cloud.points.size());
  size_t num_points = tree.size();

  for (size_t i = 0; i < num_points; ++i)
  {
    double query_pt[3];
    query_pt[0] = tree.get_x(i);
    query_pt[1] = tree.get_y(i);
    query_pt[2] = tree.get_z(i);

    size_t ret_index;
    double out_dist_sqr;

    nanoflann::KNNResultSet<double> resultSet(1);
    resultSet.init(&ret_index, &out_dist_sqr);
    index.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters());

    points_per_edge[ret_index].push_back(i);
  }

  for (size_t i = 0; i < points_per_edge.size(); ++i)
  {
    const auto& point_indices = points_per_edge[i];
    int eid = index_to_eid[i];
    auto& einfo = graph.edge(eid);
    QSMEdge& ed = einfo.data;

    if (ed.cyl_ID < 1) continue;              // cyl_ID < 1 means it is a prolongation
    if (point_indices.empty()) continue;
    if (ed.conic_allometry < 0.04) continue;  // TODO: control parameter
    if (point_indices.size() < 50) continue;

    const QSMNode& src = graph.node(einfo.source);
    const QSMNode& tgt = graph.node(einfo.target);

    utils::fitting::Vec3 axis_start = {src.x, src.y, src.z};
    utils::fitting::Vec3 axis_end   = {tgt.x, tgt.y, tgt.z};

    utils::fitting::FittingCircloid fitter;
    fitter.set_axe(axis_start, axis_end);

    for (size_t pt_idx : point_indices)
    {
      double x = tree.get_x(pt_idx);
      double y = tree.get_y(pt_idx);
      double z = tree.get_z(pt_idx);
      fitter.add_point(x, y, z);
    }

    utils::fitting::FittingResult res = fitter.fit(0.03);

    double r_meas    = res.radius;
    double p_inside  = 0;
    double arc       = res.arc_coverage_deg;
    double p_inlier  = res.inlier_percentage;
    bool valid = true;
    if (r_meas < srmeas) valid = false;
    else if (p_inside > sins) valid = false;
    else if (!(arc > sarc && p_inlier > sinl)) valid = false;

    // We fitted a valid circloid:
    //  - update the radius
    //  - update the position of the node for more accuracy
    if (valid)
    {
      ed.radius = r_meas;
      /*auto nodeid = einfo.source;
      auto& node = graph.nodes()[nodeid];
      node.x = res.center.x;
      node.y = res.center.y;
      node.z = res.center.z;*/
    }
  }
}

void QSMbuilder::refine_radii(const PointCloud& tree)
{
  ServiceLocator::logger()("Refine radii");

    if (graph.edge_count() == 0) return;

  // Prepare centroids for KD-Tree
  SimpleAdaptor centroids_cloud;
  centroids_cloud.points.reserve(graph.edge_count());
  std::vector<int> index_to_eid;
  index_to_eid.reserve(graph.edge_count());

  for (auto& [eid, einfo] : graph.edges())
  {
    const QSMNode& src = graph.node(einfo.source);
    const QSMNode& tgt = graph.node(einfo.target);
    double cx = (src.x + tgt.x) / 2.0;
    double cy = (src.y + tgt.y) / 2.0;
    double cz = (src.z + tgt.z) / 2.0;
    centroids_cloud.points.push_back({cx, cy, cz, eid});
    index_to_eid.push_back(eid);
  }

  CentroidKDTree index(3, centroids_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  // Assign tree points to nearest edge (cylinder)
  std::vector<std::vector<size_t>> points_per_edge(centroids_cloud.points.size());
  size_t num_points = tree.size();
  for (size_t i = 0; i < num_points; ++i)
  {
    double query_pt[3] = { tree.get_x(i), tree.get_y(i), tree.get_z(i) };
    size_t ret_index;
    double out_dist_sqr;
    nanoflann::KNNResultSet<double> resultSet(1);
    resultSet.init(&ret_index, &out_dist_sqr);
    index.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters());
    points_per_edge[ret_index].push_back(i);
  }

  for (size_t i = 0; i < points_per_edge.size(); ++i)
  {
    const auto& point_indices = points_per_edge[i];
    int eid = index_to_eid[i];
    auto& einfo = graph.edge(eid);
    QSMEdge& ed = einfo.data;

    // Skip too small radii. They are unlikely to be robust enough
    if (ed.radius < 0.03) continue;
    if (point_indices.empty()) continue;

    // Get the orientation of the edge by gathering its nodes
    const QSMNode& src = graph.node(einfo.source);
    const QSMNode& tgt = graph.node(einfo.target);
    utils::fitting::Vec3 axis_start = {src.x, src.y, src.z};
    utils::fitting::Vec3 axis_end   = {tgt.x, tgt.y, tgt.z};

    // Fitting
    utils::fitting::FittingCircloid fitter;
    fitter.set_axe(axis_start, axis_end);
    for (size_t pt_idx : point_indices)
    {
      fitter.add_point(tree.get_x(pt_idx), tree.get_y(pt_idx), tree.get_z(pt_idx));
    }

    utils::fitting::FittingResult res = fitter.fit(0.04);

    float ratio = (res.radius - ed.radius) / ed.radius;
    // Very strict validation: we only accept a near-perfect measurement
    bool valid =
      res.arc_coverage_deg > 320.0 &&   // Close to full closed loop
      res.inlier_percentage > 80.0 &&   // Lot of inliers
      ed.radius >= 0.04 &&              // At least 4 cm radius
      ratio > -0.1;                     // Not smaller than -10% than reference

    if (valid)
    {
      ed.radius = res.radius;
      /*auto& node = graph.nodes()[einfo.source];
      node.x = res.center.x;
      node.y = res.center.y;
      node.z = res.center.z;*/
    }
  }
}

void QSMbuilder::polynomial_fitting(double tip_radius)
{
  // Group edges by axis_ID
  std::map<int, std::vector<int>> axes;   // axis_ID -> edge IDs
  for (const auto& [eid, einfo] : graph.edges())
    axes[einfo.data.axis_ID].push_back(eid);

  for (auto& [axis_id, eids] : axes)
  {
    // Sort edges root→tip by subtree_length (descending)
    std::sort(eids.begin(), eids.end(), [this](int a, int b) {
      return graph.edge_data(a).subtree_length > graph.edge_data(b).subtree_length;
    });

    // Collect valid measurements
    std::vector<std::pair<double, double>> data_points;
    data_points.reserve(eids.size());

    for (int eid : eids)
    {
      const QSMEdge& ed = graph.edge_data(eid);
      if (ed.radius != RADIUS_UNSET && ed.subtree_length != SUBTREE_LENGTH_UNSET)
      {
        double y = ed.radius - tip_radius;
        double x = ed.subtree_length;
        data_points.push_back({x, y});
      }
    }

    if (data_points.size() <= 6) continue;

    double sum_x2 = 0, sum_x3 = 0, sum_x4 = 0, sum_xy = 0, sum_x2y = 0;

    for (const auto& p : data_points)
    {
      double x = p.first;
      double y = p.second;
      double x2 = x * x;
      double x3 = x2 * x;
      double x4 = x3 * x;

      sum_x2  += x2;
      sum_x3  += x3;
      sum_x4  += x4;
      sum_xy  += x * y;
      sum_x2y += x2 * y;
    }

    double det = sum_x2 * sum_x4 - sum_x3 * sum_x3;
    if (std::abs(det) < 1e-12) continue;

    double a = (sum_xy * sum_x4 - sum_x2y * sum_x3) / det;
    double b = (sum_x2 * sum_x2y - sum_xy * sum_x3) / det;

    std::vector<double> predictions;
    predictions.reserve(eids.size());
    for (int eid : eids)
    {
      double len = graph.edge_data(eid).subtree_length;
      double pred = tip_radius + a * len + b * len * len;
      if (pred < 0) pred = 0;
      predictions.push_back(pred);
    }

    // Ensure no increasing diameters root→tip
    double previous_radius = predictions.back();
    for (auto it = predictions.rbegin(); it != predictions.rend(); ++it)
    {
      if (*it < previous_radius) *it = previous_radius;
      previous_radius = *it;
    }

    int i = 0;
    for (int eid : eids)
    {
      QSMEdge& ed = graph.edge_data(eid);
      double pred_radius = predictions[i++];

      if (ed.subtree_length == SUBTREE_LENGTH_UNSET)
        throw std::logic_error("subtree_length unset during polynomial fitting");

      bool should_update = false;

      if (ed.radius == RADIUS_UNSET)
      {
        should_update = true;
      }
      else
      {
        double diff = std::abs(ed.radius - pred_radius);
        //double diff = 99999; // feature disabled
        if (diff > 0.25 * pred_radius)
          should_update = true;
      }

      if (should_update)
      {
        ed.radius = pred_radius;
      }
    }
  }
}

void QSMbuilder::reconstruct_missing_radii(double tip_radius)
{
  // Group edges by branch_order
  std::map<int, std::vector<int>> by_order;
  for (const auto& [eid, einfo] : graph.edges())
    by_order[einfo.data.branch_order].push_back(eid);

  for (auto& [branch_order, eids] : by_order)
  {
    if (branch_order == 1) continue;

    // Group by axis_ID
    std::unordered_map<int, std::vector<int>> axes;
    for (int eid : eids) axes[graph.edge_data(eid).axis_ID].push_back(eid);

    for (auto& [axis_id, axe_eids] : axes)
    {
      if (axe_eids.empty()) continue;

      // Sort root->tip
      std::sort(axe_eids.begin(), axe_eids.end(), [this](int a, int b) {
        return graph.edge_data(a).subtree_length > graph.edge_data(b).subtree_length;
      });

      // The first edge is the root of the branch; find its parent edge
      int first_eid  = axe_eids[0];
      int parent_cyl_id = graph.edge_data(first_eid).parent_ID;

      // Find parent edge by cyl_ID
      int parent_eid = -1;
      for (const auto& [eid, einfo] : graph.edges())
      {
        if (einfo.data.cyl_ID == parent_cyl_id) { parent_eid = eid; break; }
      }
      if (parent_eid < 0) continue;

      const QSMEdge& parent_ed = graph.edge_data(parent_eid);
      const double r0 = parent_ed.radius * 0.85;
      const double w0 = parent_ed.subtree_length;

      // Check if reconstruction is needed (any RADIUS_UNSET in axis)
      bool need_recon = false;
      for (int eid : axe_eids)
        if (graph.edge_data(eid).radius == RADIUS_UNSET) { need_recon = true; break; }

      if (need_recon)
      {
        for (int eid : axe_eids)
        {
          QSMEdge& ed = graph.edge_data(eid);
          ed.radius = conic_allometry(tip_radius, ed.subtree_length, w0, r0);
        }
      }
      /*else
      {
        double r1    = graph.edge_data(axe_eids[0]).radius;
        double ratio = r0 / r1;
        if (ratio < 1)
        {
          for (int eid : axe_eids)
            graph.edge_data(eid).radius *= ratio*0.85;
        }
      }*/
    }
  }
}

double QSMbuilder::conic_allometry(double tip_radius, double wi, double w0, double r0) const
{
  if (w0 == 0) return tip_radius;
  const double s = std::pow(wi / w0, 1.1);
  return tip_radius + s * (r0 - tip_radius);
}

void QSMbuilder::conic_allometry(double R0, double tip_radius)
{
  if (graph.edges().size() == 1)
  {
    for (auto& [eid, einfo] : graph.edges())
    {
      einfo.data.radius = R0;
    }

    return;
  }

  double w0 = 0;
  for (const auto& [eid, einfo] : graph.edges())
  {
    if (einfo.data.parent_ID == 0)
    {
      w0 = einfo.data.subtree_length;
      break;
    }
  }

  if (w0 <= 0)
  {
    // No root found or root has zero/negative subtree_length – skip allometry
    return;
  }

  for (auto& [eid, einfo] : graph.edges())
  {
    QSMEdge& ed = einfo.data;
    double wi = ed.subtree_length;
    ed.radius          = conic_allometry(tip_radius, wi, w0, R0);
    ed.conic_allometry = ed.radius;
  }
}

}

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
#include "fitting.h"

namespace arbor::qsm {

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
  distance_to_root();
  float L = 0.0f;
  for (const auto& [edge_id, edge_info] : graph.edges())
  {
    const float dist = edge_info.data.distance_to_root;

    if (dist != arbor::qsm::DISTANCE_TO_ROOT_UNSET && dist > L)
    {
      L = dist;
    }
  }

  float H = 0.0f;
  for (size_t i = 0 ; i < tree.size() ; i++)
    if (tree.get_hag(i) > H) H = tree.get_hag(i);

  H = std::max(L, H);

  // Estimation of an expected radius based on broad allometry
  double R0 = Allometry::DBH_vs_H(H) / 2.0;

  // TODO: control parameter
  // If the estimated radius is too small we don't even try to measure the tree
  if (R0 < 0.03)
  {
    conic_allometry(R0, tip_radius);

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2);
    oss << "[Small tree allometry] H = " << H
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
    if (einfo.data.axis_id == 1 && einfo.data.radius == RADIUS_UNSET)
    {
      has_na = true;
      break;
    }
  }

  if (has_na)
  {
    std::string msg = "[No valid measure] Not a single valid measure for this tree. The QSM is a pure reconstruction based on allometry";
    ServiceLocator::logger()("\033[33m" + msg + "\033[0m");
    graph.messages.push_back(msg);
    conic_allometry(R0, tip_radius);
    return;
  }


  // Check root radius
  double Rroot;
  for (auto& [eid, einfo] : graph.edges())
  {
    if (einfo.data.axis_id == 1 && einfo.data.source == 0)
    {
      Rroot = einfo.data.radius;
      break;
    }
  }

  if (Rroot > 3*R0)
  {
    std::ostringstream oss;
    oss << "[Diameter anomaly] Measured root diameter is "
        << static_cast<int>(Rroot / R0)
        << " times greater than the expected DBH ("
        << std::fixed << std::setprecision(1) << 2 * Rroot
        << " vs. "
        << std::fixed << std::setprecision(2) << 2 * R0
        << ").";

    std::string msg = oss.str();
    ServiceLocator::logger()("\033[33m" + msg + "\033[0m");
    graph.messages.push_back(msg);

    //conic_allometry(2.0 * R0, tip_radius);
    //return;
  }

  ServiceLocator::logger()("Reconstruction");
  reconstruct_missing_radii(tip_radius);
}

void QSMbuilder::measure_radii(const PointCloud& tree, float sarc, float sins, float sinl, float srmeas)
{
  auto points_by_edge = group_points_by_edge(graph, tree);
  if (points_by_edge.empty()) throw std::runtime_error("Internal error: no cylinders found.");

  for (auto& [eid, point_indices] : points_by_edge)
  {
    auto& einfo = graph.edge(eid);
    QSMEdge& ed = einfo.data;
    ed.radius = RADIUS_UNSET;

    // Filtration logic
    if (point_indices.size() < 50 || ed.conic_allometry < 0.03)
      continue;

    const auto& src = graph.node(einfo.source);
    const auto& tgt = graph.node(einfo.target);

    utils::fitting::FittingCircloid fitter;
    fitter.set_axe({src.x, src.y, src.z}, {tgt.x, tgt.y, tgt.z});

    for (size_t pt_idx : point_indices)
      fitter.add_point(tree.get_x(pt_idx), tree.get_y(pt_idx), tree.get_z(pt_idx));

    auto res = fitter.fit(0.03);

    bool valid = (res.radius >= srmeas) && (res.arc_coverage_deg > sarc) && (res.inlier_percentage > sinl);

    if (valid)
    {
      ed.radius = res.radius;
      ed.quality = MEASURED;
    }
  }
}

void QSMbuilder::refine_radii(const PointCloud& tree)
{
  ServiceLocator::logger()("Refine radii");

  auto points_by_edge = group_points_by_edge(graph, tree);

  for (auto& [eid, point_indices] : points_by_edge)
  {
    auto& einfo = graph.edge(eid);
    if (einfo.data.radius < 0.03 || point_indices.empty()) continue;

    const auto& src = graph.node(einfo.source);
    const auto& tgt = graph.node(einfo.target);

    utils::fitting::FittingCircloid fitter;
    fitter.set_axe({src.x, src.y, src.z}, {tgt.x, tgt.y, tgt.z});

    for (size_t pt_idx : point_indices)
      fitter.add_point(tree.get_x(pt_idx), tree.get_y(pt_idx), tree.get_z(pt_idx));

    auto res = fitter.fit(0.04);
    float ratio = (res.radius - einfo.data.radius) / einfo.data.radius;

    bool valid = res.arc_coverage_deg > 300.0 && res.inlier_percentage > 70.0 && einfo.data.radius >= 0.04 && ratio > -0.1;

    if (valid)
    {
      einfo.data.radius = res.radius;
      einfo.data.quality = REFINED;
    }
  }
}

void QSMbuilder::polynomial_fitting(double tip_radius)
{
  // Group edges by axis_id
  std::map<int, std::vector<int>> axes;   // axis_id -> edge IDs
  for (const auto& [eid, einfo] : graph.edges())
    axes[einfo.data.axis_id].push_back(eid);

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
        ed.quality = POLYNOMIAL;
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

    // Group by axis_id
    std::unordered_map<int, std::vector<int>> axes;
    for (int eid : eids) axes[graph.edge_data(eid).axis_id].push_back(eid);

    for (auto& [axis_id, axe_eids] : axes)
    {
      if (axe_eids.empty()) continue;

      // Sort root->tip
      std::sort(axe_eids.begin(), axe_eids.end(), [this](int a, int b) {
        return graph.edge_data(a).subtree_length > graph.edge_data(b).subtree_length;
      });

      // The first edge is the root of the branch; find its parent edge
      int first_eid = axe_eids[0];
      int parent_id = graph.edge_data(first_eid).source;

      // Find parent edge by id
      int parent_eid = -1;
      for (const auto& [eid, einfo] : graph.edges())
      {
        if (einfo.data.id == parent_id) { parent_eid = eid; break; }
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
          ed.quality = CONICALLOM;
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
      einfo.data.quality = CONICALLOM;
    }

    return;
  }

  double w0 = 0;
  for (const auto& [eid, einfo] : graph.edges())
  {
    if (einfo.data.source == 0)
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
    ed.quality = CONICALLOM;
  }
}

}

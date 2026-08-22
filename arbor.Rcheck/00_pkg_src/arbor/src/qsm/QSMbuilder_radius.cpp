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
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>

#include "QSMbuilder.h"
#include "PointCloud.h"

#include "allometry.h"
#include "fitting.h"
#include "fitting_quality.h"
#include "fitting_polynomial.h"

namespace arbor::qsm {

bool QSMbuilder::construct_radii(const PointCloud& tree, double tip_radius)
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
  auto model = AllometryDataBase::getAllometry(params.qsm.allometry_name);
  float D0 = model->DBH_vs_H(H);
  D0 *= params.qsm.allometry_scale;
  float R0 = D0/2.0f;

  // If the estimated radius is too small we don't even try to measure the tree
  if (D0 < params.qsm.min_measurable_dbh && !likely_broken)
  {
    conic_allometry(R0, tip_radius);

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2);
    oss << "[Small tree allometry] H = " << H
        << std::fixed << std::setprecision(1)
        << " m: estimated DBH = " << (D0 * 100)
        << " cm. Too small to be measured. "
        << "The QSM is the result of pure allometry without actual measurements.";

    std::string msg = oss.str();

    graph.messages.push_back(msg);
    ServiceLocator::logger()("\033[33m" + msg + "\033[0m");
    return true;
  }

  if (D0 < params.qsm.min_measurable_dbh && likely_broken)
  {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2);
    oss << "[Broken tree] H = " << H
        << std::fixed << std::setprecision(1)
        << " m: estimated DBH = " << (D0 * 100)
        << " cm. Too small to be measured but with good measurement anyway. "
        << "The QSM is likely a dead tree.";

    std::string msg = oss.str();

    graph.messages.push_back(msg);
    ServiceLocator::logger()("\033[33m" + msg + "\033[0m");
  }
  else
  {
    likely_broken = false;
  }

  ServiceLocator::logger()("Pre-allometry");
  conic_allometry(1.5*R0, tip_radius); // Overestimate the tree on purpose

  ServiceLocator::logger()("Measuring diameters");
  measure_radii(tree);

  for (auto& e : graph.edges())
  {
    if (e.second.data.quality == CONICALLOM)
    {
      e.second.data.radius  = RADIUS_UNSET;
      e.second.data.quality = UNKNOWN;
    }
  }

  compute_architecture(true);

  ServiceLocator::logger()("Polynomial fitting");
  polynomial_fitting(tip_radius);
  polynomial_fitting_root();

  bool all_na = true;
  for (auto& [eid, einfo] : graph.edges())
  {
    if (einfo.data.axis_id == 1 && einfo.data.radius != RADIUS_UNSET)
    {
      all_na = false;
      break;
    }
  }

  if (all_na)
  {
    conic_allometry(R0, tip_radius);
    return false;
  }

  pipe_model_reconstruction(tip_radius);

  // Check if main axis has valid measurements
  /*bool has_na = false;
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
    conic_allometry(R0, tip_radius);
    return false;
  }*/


  // Check root radius
  double Rroot = 0;
  for (auto& [eid, einfo] : graph.edges())
  {
    if (einfo.data.axis_id == 1 && einfo.data.source == 0)
    {
      Rroot = einfo.data.radius;
      break;
    }
  }

  if (Rroot > 3*R0 && !likely_broken)
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
  //reconstruct_missing_radii(tip_radius);

  return true;
}

void QSMbuilder::measure_radii(const PointCloud& tree)
{
  auto points_by_edge = group_points_by_edge(tree);
  if (points_by_edge.empty()) throw std::runtime_error("Internal error: no cylinders found.");

  for (auto& [eid, point_indices] : points_by_edge)
  {
    // This happens if a point is not assigned an edge id (too far?)
    if (eid <= 0) continue;

    // This happens if some edges were removed in:
    // - prune_spurious_branches()
    // - clean_tree_butt()
    // - bug in the code but it becomes invisible
    if (!graph.has_edge(eid)) continue;

    auto& einfo = graph.edge(eid);
    QSMEdge& ed = einfo.data;
    ed.radius = RADIUS_UNSET;

    // Filtration logic
    if (point_indices.size() < 50 || ed.conic_allometry < params.qsm.min_measurable_radius)
      continue;

    const auto& src = graph.node(einfo.source);
    const auto& tgt = graph.node(einfo.target);

    utils::fitting::FitQuality fit_quality  = utils::fitting::FitQuality::standard_preset();
    fit_quality.min_radius = params.qsm.min_measurable_radius;
    utils::fitting::FitQuality fit_hquality = utils::fitting::FitQuality::accurate_preset();

    utils::fitting::CrossSectionFitter fitter;
    fitter.set_axis({src.x, src.y, src.z}, {tgt.x, tgt.y, tgt.z});

    for (size_t pt_idx : point_indices)
      fitter.add_point(tree.get_x(pt_idx), tree.get_y(pt_idx), tree.get_z(pt_idx));

    utils::fitting::FitMode strategy = utils::fitting::FitMode::Standard;
    if (ed.conic_allometry > 0.15) strategy = utils::fitting::FitMode::Buttress;

    auto res = fitter.fit(fit_quality.ransac_tolerance, strategy);

    if (fit_quality.accept(res))
    {
      ed.radius = res.radius;

      if (fit_hquality.accept(res))
        ed.quality = REFINED;
      else
        ed.quality = MEASURED;
    }
  }
}

void QSMbuilder::refine_radii_broken(const PointCloud& tree)
{
  ServiceLocator::logger()("Refine radii (broken tree)");

  auto points_by_edge = group_points_by_edge(tree);

  // Group edges by axis_id, then process each axis root->tip
  std::map<int, std::vector<int>> axes;
  for (const auto& [eid, einfo] : graph.edges())
    axes[einfo.data.axis_id].push_back(eid);

  constexpr utils::fitting::FitQuality fit_quality = utils::fitting::FitQuality::accurate_and_large_preset();

  for (auto& [axis_id, eids] : axes)
  {
    // Sort root->tip by decreasing subtree_length
    std::sort(eids.begin(), eids.end(), [this](int a, int b)
    {
      return graph.edge_data(a).subtree_length > graph.edge_data(b).subtree_length;
    });

    // Stop processing an axis once 5 consecutive edges fail to reach MEASURED quality
    int consecutive_below_measured = 0;

    for (int eid : eids)
    {
      if (consecutive_below_measured >= 5) break;

      auto& einfo = graph.edge(eid);
      auto it = points_by_edge.find(eid);

      if (it == points_by_edge.end() || it->second.empty())
      {
        if (einfo.data.quality < MEASURED)
          consecutive_below_measured++;
        else
          consecutive_below_measured = 0;

        continue;
      }

      const auto& point_indices = it->second;
      const auto& src = graph.node(einfo.source);
      const auto& tgt = graph.node(einfo.target);

      utils::fitting::CrossSectionFitter fitter;
      fitter.set_axis({src.x, src.y, src.z}, {tgt.x, tgt.y, tgt.z});

      for (size_t pt_idx : point_indices)
        fitter.add_point(tree.get_x(pt_idx), tree.get_y(pt_idx), tree.get_z(pt_idx));

      auto res = fitter.fit(fit_quality.ransac_tolerance, utils::fitting::FitMode::Standard);

      const float current_radius = einfo.data.radius;

      if (fit_quality.accept(res, current_radius))
      {
        einfo.data.radius  = res.radius;
        einfo.data.quality = REFINED;
        consecutive_below_measured = 0;
      }
      else
      {
        if (einfo.data.quality < MEASURED)
          consecutive_below_measured++;
        else
          consecutive_below_measured = 0;
      }
    }
  }
}

void QSMbuilder::refine_radii(const PointCloud& tree)
{
  ServiceLocator::logger()("Refine radii");

  auto points_by_edge = group_points_by_edge(tree);

  constexpr utils::fitting::FitQuality fit_quality = utils::fitting::FitQuality::accurate_preset();

  for (auto& [eid, point_indices] : points_by_edge)
  {
    // This happens if a point is not assigned an edge id (too far?)
    if (eid <= 0) continue;

    // This happens if some edges were removed in:
    // - prune_spurious_branches()
    // - clean_tree_butt()
    // - bug in the code but it becomes invisible
    if (!graph.has_edge(eid)) continue;

    auto& einfo = graph.edge(eid);
    if (einfo.data.radius < params.qsm.min_measurable_radius || point_indices.empty()) continue;

    const auto& src = graph.node(einfo.source);
    const auto& tgt = graph.node(einfo.target);

    utils::fitting::CrossSectionFitter fitter;
    fitter.set_axis({src.x, src.y, src.z}, {tgt.x, tgt.y, tgt.z});

    for (size_t pt_idx : point_indices)
      fitter.add_point(tree.get_x(pt_idx), tree.get_y(pt_idx), tree.get_z(pt_idx));

    auto res = fitter.fit(fit_quality.ransac_tolerance, utils::fitting::FitMode::Standard);

    if (fit_quality.accept(res, einfo.data.radius))
    {
      einfo.data.radius = res.radius;
      einfo.data.quality = REFINED;
    }
  }
}

void QSMbuilder::polynomial_fitting(double tip_radius)
{
  auto axes = QSMbuilder::build_axis_map();

  for (auto& [axis_id, eids] : axes)
  {
    std::vector<std::pair<float, float>> data_points;
    data_points.reserve(eids.size());
    for (int eid : eids)
    {
      const QSMEdge& ed = graph.edge_data(eid);
      if (ed.radius != RADIUS_UNSET && ed.subtree_length != SUBTREE_LENGTH_UNSET && ed.distance_to_root > 1.3)
        data_points.push_back({ed.subtree_length, ed.radius - tip_radius});
    }

    // Don't fit a polynomial with less than 6 radius measures
    if (data_points.size() <= 6)
    {
      // Erase previous data
      for (int eid : eids)
      {
        QSMEdge& ed = graph.edge_data(eid);
        ed.radius = RADIUS_UNSET;
        ed.quality = UNKNOWN;
      }
      continue;
    }

    // The last point in data_points represents the latest valid measurement limit.
    // data_points.back().first holds the subtree_length of this boundary.
    float latest_valid_subtree_length = data_points.back().first;

    auto mode = fitting::PolynomialFitting::Mode::DecreasingOnly;
    fitting::PolynomialFitting fit(data_points, tip_radius, mode);
    if (!fit.valid) continue;

    for (int eid : eids)
    {
      QSMEdge& ed = graph.edge_data(eid);
      if (ed.subtree_length == SUBTREE_LENGTH_UNSET)
        throw std::logic_error("subtree_length unset during polynomial fitting");

      // Check if we are within the limit of the latest valid measurement
      // (Depending on your coordinate system, this might be <= or >=.
      // Assuming subtree_length decreases towards the tips, adjust the comparison if needed.)
      if (ed.subtree_length >= latest_valid_subtree_length)
      {
        double pred_radius = fit.predict(ed.subtree_length);

        // Replace previous measure only if
        // 1. There was no measure (interpolation and extension)
        // 2. There was a measure that is significantly different (30%) and NOT of high quality
        if ((ed.radius == RADIUS_UNSET) ||
            ((std::abs(ed.radius - pred_radius) > 0.30 * pred_radius) && (ed.quality < REFINED)))// && ed.distance_to_root > 1.3))
        {
          ed.radius  = pred_radius;
          ed.quality = POLYNOMIAL;
        }
      }
    }
  }
}

void QSMbuilder::polynomial_fitting_root()
{
  auto axes = QSMbuilder::build_axis_map();

  for (auto& [axis_id, eids] : axes)
  {
    if (axis_id != 1) continue;

    std::vector<std::pair<float, float>> data_points;
    data_points.reserve(eids.size());
    for (int eid : eids)
    {
      const QSMEdge& ed = graph.edge_data(eid);
      if (ed.radius != RADIUS_UNSET && ed.subtree_length != SUBTREE_LENGTH_UNSET && ed.distance_to_root <= 1.3)
        data_points.push_back({ed.subtree_length, ed.radius});
    }

    // Don't fit a polynomial with less than 4 radius measures
    if (data_points.size() <= 4) continue;

    // The last point in data_points represents the latest valid measurement limit.
    // data_points.back().first holds the subtree_length of this boundary.
    float latest_valid_subtree_length = data_points.back().first;

    auto mode = fitting::PolynomialFitting::Mode::DecreasingOnly;
    fitting::PolynomialFitting fit(data_points, 0, mode);
    if (!fit.valid) continue;

    for (int eid : eids)
    {
      QSMEdge& ed = graph.edge_data(eid);
      if (ed.subtree_length == SUBTREE_LENGTH_UNSET)
        throw std::logic_error("subtree_length unset during polynomial fitting");

      // Check if we are within the limit of the latest valid measurement
      if (ed.subtree_length >= latest_valid_subtree_length)
      {
        double pred_radius = fit.predict(ed.subtree_length);

        if ((ed.radius == RADIUS_UNSET) ||
            ((std::abs(ed.radius - pred_radius) > 0.30 * pred_radius) && (ed.quality < REFINED) && ed.distance_to_root <= 1.3))
        {
          ed.radius  = pred_radius;
          ed.quality = POLYNOMIAL;
        }
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
      const double r0 = parent_ed.radius * 0.75;
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

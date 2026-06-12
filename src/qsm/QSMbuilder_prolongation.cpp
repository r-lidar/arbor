/**
 * @file QSMbuilder_prolongation.cpp
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

#include "QSMbuilder.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace arbor::qsm {

void QSMbuilder::prolongate(double d, double L)
{
  if (d <= 0.0) return;

  ServiceLocator::logger()("Prolongation to the ground: " + std::to_string(d) + " m");

  // Collect main axis edges (axis_id == 1), ordered root → tip (descending subtree_length)
  std::vector<int> axis_eids;
  for (const auto& [eid, einfo] : graph.edges())
  {
    if (einfo.data.axis_id == 1)
      axis_eids.push_back(eid);
  }

  const size_t n = axis_eids.size();
  if (n < 2) return;

  std::sort(axis_eids.begin(), axis_eids.end(), [this](int a, int b) {
    return graph.edge_data(a).subtree_length > graph.edge_data(b).subtree_length;
  });

  int root_eid     = axis_eids[0];
  double root_radius = graph.edge_data(root_eid).radius;

  // Compute cylinder lengths and cumulative lengths
  std::vector<double> lens(n);
  std::vector<double> cum(n);
  double total = 0.0;
  for (size_t i = 0; i < n; ++i)
  {
    const auto& einfo = graph.edge(axis_eids[i]);
    lens[i] = graph.edge_data(axis_eids[i]).length(
      graph.node(einfo.source), graph.node(einfo.target));
    total += lens[i];
    cum[i] = total;
  }
  if (total == 0.0) return;

  // Take only the first min(3.0 m, 25% of total axis length) for direction estimation
  double cutoff = std::min(3.0, 0.25 * total);
  size_t k = 0;
  while (k < n && cum[k] <= cutoff) k++;
  if (k == 0) k = 1;
  if (k > n)  k = n;

  NodeID root_src_nid = graph.edge(root_eid).source;
  const QSMNode& root_node = graph.node(root_src_nid);
  const QSMNode& last_node = graph.node(graph.edge(axis_eids[k - 1]).target);

  // Estimated orientation of the trunk (root -> tip direction)
  double dx = last_node.x - root_node.x;
  double dy = last_node.y - root_node.y;
  double dz = last_node.z - root_node.z;
  double N  = std::sqrt(dx*dx + dy*dy + dz*dz);
  if (N <= 0.0) return;

  double ox = dx / N;
  double oy = dy / N;
  double oz = dz / N;

  double angle_deg = std::acos(std::clamp(oz, -1.0, 1.0)) * 180.0 / M_PI;

  if (angle_deg > 50.0)
  {
    constexpr double MAX_ANGLE_RAD = 50.0 * M_PI / 180.0;
    double horiz = std::sqrt(ox*ox + oy*oy);  // current horizontal magnitude

    oz = std::cos(MAX_ANGLE_RAD);             // cos(50°) ≈ 0.643
    double horiz_clamped = std::sin(MAX_ANGLE_RAD); // sin(50°) ≈ 0.766

    ox = (ox / horiz) * horiz_clamped;
    oy = (oy / horiz) * horiz_clamped;

    auto angle_str = std::to_string(angle_deg);
    angle_str.resize(angle_str.find('.') + 3); // 2 decimals
    ServiceLocator::logger()("Trunk angle " + angle_str + " deg exceeds 50 deg, clamping to 50 deg for prolongation");
  }

  double d_adj;
  if (oz < 1e-9)
    d_adj = d;
  else
    d_adj = d / oz;

  char buf[256];
  std::snprintf(buf, sizeof(buf), "Trunk angle: %.2f deg | Distance %.2f m adjusted to %.2f m", angle_deg, d, d_adj);
  ServiceLocator::logger()(buf);

  int nseg = std::max(1, int(std::ceil(d_adj / L)));

  double root_subtree = graph.edge_data(root_eid).subtree_length;
  if (root_subtree == SUBTREE_LENGTH_UNSET)
    throw std::runtime_error("Invalid QSM, the root has no subtree length");

  constexpr int NO_PARENT = -1;

  int next_id = -nseg;          // IDs run -nseg … -1 for the new segments
  int prev_id = NO_PARENT;      // The deepest segment has no parent
  NodeID prev_node_id = -1;

  // Build from deepest (new root, i = nseg) up to shallowest (i = 1)
  for (int i = nseg; i >= 1; i--)
  {
    double f1 = double(i)     / nseg;   // fraction for the deeper endpoint
    double f2 = double(i - 1) / nseg;   // fraction for the shallower endpoint

    double x1 = root_node.x - ox * d_adj * f1;
    double y1 = root_node.y - oy * d_adj * f1;
    double z1 = root_node.z - oz * d_adj * f1;

    double x2 = root_node.x - ox * d_adj * f2;
    double y2 = root_node.y - oy * d_adj * f2;
    double z2 = root_node.z - oz * d_adj * f2;

    NodeID node_id_1, node_id_2;

    if (i == nseg)
    {
      // Deepest node: create the new root node
      node_id_1 = graph.add_node({x1, y1, z1});
    }
    else
    {
      // Reuse the previously created shallower node as this segment's deep end
      node_id_1 = prev_node_id;
    }

    if (i == 1)
    {
      // Shallowest segment: connect directly to the original root node
      node_id_2 = root_src_nid;
    }
    else
    {
      node_id_2 = graph.add_node({x2, y2, z2});
    }

    QSMEdge ed;
    ed.id           = next_id;
    ed.source       = prev_id;
    ed.axis_id      = 1;
    ed.branch_order = 1;
    ed.subtree_length = root_subtree + (d_adj / nseg) * i;
    ed.radius = root_radius * std::pow(1.01, i - 1);

    graph.add_edge(node_id_1, node_id_2, ed);

    prev_id      = next_id;
    prev_node_id = node_id_2;
    next_id++;
  }

  // Connect the original root edge to the prolongation chain.
  // prev_id is now -1, the ID of the shallowest new segment.
  graph.edge_data(root_eid).source = prev_id;
}

void QSMbuilder::estimate_prolongation(const PointCloud& tree)
{
  prolongation_distance = 0.0;

  if (!tree.has_hag() || graph.edge_count() == 0) return;

  // Compute the minimum startZ across all edges (source node Z) –
  // matches the original behaviour that used min(startZ) over all cylinders.
  double min_start_z = std::numeric_limits<double>::max();
  for (const auto& [eid, einfo] : graph.edges())
  {
    double sz = graph.node(einfo.source).z;
    if (sz < min_start_z) min_start_z = sz;
  }

  // Find max(hag) for points where Z <= min_start_z
  double max_hag = -std::numeric_limits<double>::max();
  bool found = false;

  for (size_t i = 0; i < tree.size(); ++i)
  {
    if (tree.get_z(i) <= min_start_z)
    {
      double hag = tree.get_hag(i);
      if (hag > max_hag) { max_hag = hag; found = true; }
    }
  }

  if (found) prolongation_distance = max_hag;
}

}

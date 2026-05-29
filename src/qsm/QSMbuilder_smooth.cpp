/**
 * @file QSMbuilder_smooth.cpp
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
#include <map>

namespace arbor::qsm {

// Apply Taubin smoothing in-place on a 1-D signal.
// Endpoints are fixed (index 0 and index n-1 are never moved).
// lambda > 0 : shrink pass
// mu     < 0 : inflate pass  (|mu| > lambda biases toward inflation)
static void taubin_smooth(std::vector<double>& v, int steps, double lambda, double mu,  double floor_value = 0.0)
{
  const int n = static_cast<int>(v.size());
  if (n < 3 || steps == 0) return;

  for (int s = 0; s < steps; ++s)
  {
    // Pass 1: Shrink
    for (int i = 1; i < n - 1; ++i)
    {
      double delta = 0.5 * (v[i-1] + v[i+1]) - v[i];
      v[i] += lambda * delta;
      if (floor_value > 0.0) v[i] = std::max(v[i], floor_value);
    }

    // Pass 2: Inflate
    for (int i = 1; i < n - 1; ++i)
    {
      double delta = 0.5 * (v[i-1] + v[i+1]) - v[i];
      v[i] += mu * delta;
      if (floor_value > 0.0) v[i] = std::max(v[i], floor_value);
    }
  }
}

void QSMbuilder::smooth_skeleton(int steps)
{
  if (steps == 0 || graph.edge_count() == 0) return;

  ServiceLocator::logger()("Smoothing skeleton");

  for (const auto& [axis_id, edge_ids] : build_axis_map())
  {
    const int n = static_cast<int>(edge_ids.size());
    if (n < 3) continue;  // Need at least 3 edges to smooth

    // 'steps' is reduced by 3/4 at each branch order to avoid
    // over-smoothing small-diameter branches (no diameters yet at this stage).
    int branch_order = graph.edge(edge_ids[0]).data.branch_order;
    int curr_steps = static_cast<int>(std::floor(steps * std::pow(3.0/4.0, static_cast<double>(branch_order - 1))));
    curr_steps = std::max(curr_steps, 2);

    // Collect node coordinates along the axis (source nodes + last target).
    std::vector<double> x(n + 1), y(n + 1), z(n + 1);
    for (int i = 0; i < n; ++i)
    {
      const QSMNode& node = graph.node(graph.edge(edge_ids[i]).source);
      x[i] = node.x;  y[i] = node.y;  z[i] = node.z;
    }
    const QSMNode& tip = graph.node(graph.edge(edge_ids[n - 1]).target);
    x[n] = tip.x;  y[n] = tip.y;  z[n] = tip.z;

    taubin_smooth(x, curr_steps, 0.5, -0.53);
    taubin_smooth(y, curr_steps, 0.5, -0.53);
    taubin_smooth(z, curr_steps, 0.5, -0.53);

    // Write smoothed coordinates back.
    for (int i = 0; i < n; ++i)
    {
      QSMNode& node = graph.node(graph.edge(edge_ids[i]).source);
      node.x = x[i];  node.y = y[i];  node.z = z[i];
    }
    QSMNode& tip_ = graph.node(graph.edge(edge_ids[n - 1]).target);
    tip_.x = x[n];  tip_.y = y[n];  tip_.z = z[n];
  }
}

void QSMbuilder::smooth_radii()
{
  if (graph.edge_count() == 0) return;

  ServiceLocator::logger()("Smoothing radii");

  for (const auto& [axis_id, edge_ids] : build_axis_map())
  {
    const int n = static_cast<int>(edge_ids.size());
    if (n < 3) continue;

    std::vector<double> r(n);
    for (int i = 0; i < n; ++i)
      r[i] = graph.edge_data(edge_ids[i]).radius;

    // Floor: half the smallest radius on the axis, with an absolute minimum.
    double r_min = *std::min_element(r.begin(), r.end()) * 0.5;
    r_min = std::max(r_min, 1e-4);

    // Root radius (index 0) is fixed inside taubin_smooth.
    // Tip (index n-1) is also fixed, so restore it after smoothing.
    double r_tip = r[n - 1];

    taubin_smooth(r, 15, 0.5, -0.7, r_min);

    r[n - 1] = r_tip;  // Re-pin the tip radius

    for (int i = 0; i < n; ++i)
      graph.edge(edge_ids[i]).data.radius = r[i];
  }
}

}

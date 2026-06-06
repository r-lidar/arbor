/**
 * @file QSMbuilder.cpp
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

#include <limits>
#include <vector>
#include <unordered_set>

namespace arbor::qsm {

static std::vector<int> cut(const std::vector<float>& x, float by = 0.1);

void QSMbuilder::build(const PointCloud& tree)
{
  std::size_t n = tree.size();

  if (n == 0) throw std::runtime_error("Empty point cloud");

  // Find the vertical bounds of the tree
  double min_z = std::numeric_limits<double>::max();
  double max_z = -std::numeric_limits<double>::max();

  std::size_t index;
  for (std::size_t i = 0; i < n; ++i)
  {
    double current_z = tree.get_z(i);
    if (current_z < min_z) min_z = current_z;
    if (current_z > max_z)
    {
      max_z = current_z;
      index = i;
    }
  }

  // Calculate height and 1% threshold
  double tree_height = max_z - min_z;
  double trim_offset = tree_height * 0.01;
  double z_threshold = min_z + trim_offset;

  // Filter wood only
  // Cut below the 1% threshold to clean the bottom of the tree
  int count = 0;
  std::vector<bool> wood_mask(n, false);
  for (std::size_t i = 0; i < n; ++i)
  {
    bool b = tree.is_wood(i) && (tree.get_z(i) >= z_threshold);
    if (b) count++;
    wood_mask[i] = b;
    if (i == index) wood_mask[i] = true; // The highest point is enforced to be wood.
  }

  if (count == 0)
  {
    graph.messages.push_back("[No wood point] This tree has no point labelled as wood");
    return;
  }

  PointCloud wood = tree.subset(wood_mask);
  n = wood.size();

  if (n == 0) throw std::runtime_error("No wood points");

  // Determine the geographic coordinates minimum
  // and center on 0,0,0 for numerical stability
  size_t min_idx = 0;
  min_z = wood.get_z(0);
  for (size_t i = 0; i < n; ++i)
  {
    double current_z = wood.get_z(i);
    if (current_z < min_z)
    {
      min_z = current_z;
      min_idx = i;
    }
  }
  double tx = wood.get_x(min_idx);
  double ty = wood.get_y(min_idx);
  double tz = wood.get_z(min_idx);
  wood.translate(tx, ty, tz);

  // Sometime the very bottom have two clusters of wood
  // this creates troubles. One of the two is removed to ensure
  // a single entry point for the QSM
  wood = clean_tree_butt(wood);
  n = wood.size();

  // We need ground points. We don't have ground points so
  // as a workaround we use first 5 cm above min HAG
  PointCloud gnd;
  double th = 0.05;
  while (gnd.size() == 0)
  {
    std::vector<bool> gnd_mask(n);
    for (size_t i = 0; i < n; i++)
    {
      if (wood.get_z(i) < th)
        gnd_mask[i] = true;
    }
    gnd = wood.subset(gnd_mask);
    th += 0.05;
  }

  if (gnd.size() == wood.size())
    throw std::runtime_error("Internal error in QSMbuilder::build: gnd size == wood size. Please report.");

  arbor::settings::GraphParameters p;
  p.k = 80;
  p.max_gap = 1.5;
  p.power = 1.5;
  std::vector<float> dist = arbor::segment::dist2root(wood, gnd, p);

  // Keep only points with valid distance to root
  std::vector<bool> valid_mask(n);
  size_t n_valid = 0;
  for (size_t i = 0; i < n; ++i)
  {
    if (dist[i] >= 0)
    {
      valid_mask[i] = true;
      ++n_valid;
    }
  }

  wood = wood.subset(valid_mask);

  std::vector<float> filtered_dist;
  filtered_dist.reserve(n_valid);

  for (size_t i = 0; i < n; ++i)
  {
    if (valid_mask[i])
      filtered_dist.push_back(dist[i]);
  }
  dist = std::move(filtered_dist);
  n = wood.size();

  // Now we can compute iter and clust.
  // What is iter and cluster varies depending on arbor versions
  std::vector<int>    iter = cut(dist, params.qsm.skeleton_node_distance);
  std::vector<int>   clust = cluster(wood, iter, params.qsm.dbscan_eps_distance);

  std::vector<std::pair<int, int>> iter_cluster;
  iter_cluster.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    iter_cluster.emplace_back(iter[i], clust[i]);
  }

  // Build the QSM nodes
  auto ids = build_skeleton(wood, iter_cluster, params.qsm.max_d);

  // Assign cluster IDs to the point cloud (used later)
  for (size_t i = 0 ; i < wood.size() ; i ++)
    wood.set_treeid(i, ids[i]);

  ids.clear();
  ids.shrink_to_fit();

  // ============================

  // Extremely rare case with so few points that we have no cluster
  // (seen once with a very bad DTM in Murray's data)
  if (graph.edges().size() == 0)
  {
    shift(tx, ty, tz);
    return;
  }

  compute_topology();

  // Fix root issue (rare)
  int n_root = count_nodes_connected_to_root();
  if (n_root == 0) throw std::runtime_error("Internal error in QSMbuilder::build. 0 root for this QSM. Please report.");
  if (n_root > 1)
  {
    ServiceLocator::logger()("Multiple nodes connected to root detected");
    fix_multiple_root();
  }

  compute_architecture(false);
  smooth_skeleton(params.qsm.smooth_steps);

  if (likely_broken)  ServiceLocator::logger()("Detection of likely broken tree");

  detect_weird_butt();
  construct_radii(wood, params.qsm.apex_radius);
  compute_architecture(true);
  smooth_skeleton(params.qsm.smooth_steps/2);

  if (!likely_broken || !params.qsm.broken_detection_enabled)
  {
    // Deleted the step. contruct radii now flag the quality of fit
    //refine_radii(wood);
  }
  else
  {
    refine_radii_broken(wood);
    compute_architecture(true);
    smooth_skeleton(params.qsm.smooth_steps/2);
  }

  estimate_prolongation(wood);
  prolongate(prolongation_distance);
  smooth_radii();
  distance_to_root();
  shift(tx, ty, tz);

  graph.validate();
}

// Build axis_id -> edge IDs map, with each vector sorted root-to-tip (descending subtree_length).
std::map<int, std::vector<int>>  QSMbuilder::build_axis_map()
{
  std::map<int, std::vector<int>> axis_map;
  for (const auto& [eid, einfo] : graph.edges())
    axis_map[einfo.data.axis_id].push_back(eid);

  for (auto& [axis, vec] : axis_map)
  {
    std::sort(vec.begin(), vec.end(), [this](int a, int b)
    {
      return graph.edge_data(a).subtree_length > graph.edge_data(b).subtree_length;
    });
  }

  return axis_map;
}

std::unordered_map<int, std::vector<size_t>> QSMbuilder::group_points_by_edge(const PointCloud& tree)
{
  std::unordered_map<int, std::vector<size_t>> points_per_eid;

  for (size_t i = 0 ; i < tree.size() ; i++)
    points_per_eid[tree.get_treeid(i)].push_back(i);

  return points_per_eid;
}

std::vector<int> cut(const std::vector<float>& x, float by)
{
  if (x.empty()) return {};

  float xmax = *std::max_element(x.begin(), x.end());
  float xmin = *std::min_element(x.begin(), x.end());

  float start = std::floor(xmin);
  float end   = std::ceil(xmax);

  const float eps = 1e-9;

  int n_bins = static_cast<int>(std::round((end - start) / by));

  std::vector<int> out;
  out.reserve(x.size());

  for (double val : x)
  {
    int bin = static_cast<int>(std::floor(((val - start) + eps) / by));
    if (bin >= n_bins) bin = n_bins - 1;
    if (bin < 0) bin = 0;
    out.push_back(bin);
  }

  return out;
}

}


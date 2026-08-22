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
static std::pair<int, int> which_min_max_z(const PointCloud& tree);
static PointCloud make_ground_seed(const PointCloud&);
static PointCloud wood_subset(const PointCloud&, double z_threshold, int max_z_idx);

void QSMbuilder::build(const PointCloud& tree)
{
  bool retried = false;

  std::size_t n = tree.size();

  if (n == 0) throw std::runtime_error("Empty point cloud");

  // Find tree bounds
  auto min_max_idx = which_min_max_z(tree);
  double min_z = tree.get_z(min_max_idx.first);
  double max_z = tree.get_z(min_max_idx.second);

  // Calculate height and 1% threshold
  double tree_height = max_z - min_z;
  double trim_offset = tree_height * 0.01;
  double z_threshold = min_z + trim_offset;

  // Filter wood only +
  // Cut below the 1% threshold to clean the bottom of the tree +
  // Force to include the highest point as part of wood
  PointCloud wood = wood_subset(tree, z_threshold, min_max_idx.second);
  n = wood.size();

  if (n == 0) {
    graph.messages.push_back("[No wood point] This tree has no point labelled as wood");
    return;
  }

  // Determine the geographic coordinates minimum
  // and center on 0,0,0 for numerical stability
  min_max_idx = which_min_max_z(wood);
  int    min_idx = min_max_idx.first;
  double tx = wood.get_x(min_idx);
  double ty = wood.get_y(min_idx);
  double tz = wood.get_z(min_idx);
  wood.translate(tx, ty, tz);

  // Sometime the very bottom have two clusters of wood
  // this creates troubles. One of the two is removed to ensure
  // a single entry point for the QSM
  wood = clean_tree_butt(wood);
  n = wood.size();

  // We need ground points to resolve the graph.
  // We don't have ground points so, as a workaround,
  // we use first 5 cm above 0. The tree is already centered
  // on (0,0,0)
  PointCloud gnd = make_ground_seed(wood);

  // Resolve the graph to get the distance to the ground for each point
  // The distance to the ground is how we make layers in the tree.
  arbor::settings::GraphParameters p;
  p.k       = 80;
  p.max_gap = 1.5;
  p.power   = 1.5;
  std::vector<float> dist = arbor::segment::dist2root(wood, gnd, p);

  // Keep only points with valid distance to root. Some may have a distance = -1
  // meaning they were not reachable (too far from other points)
  std::vector<bool> valid_mask(n);
  std::transform(dist.begin(), dist.end(), valid_mask.begin(), [](float d) { return d >= 0.f; });
  wood = wood.subset(valid_mask);
  n = wood.size();

  std::vector<float> filtered_dist;
  filtered_dist.reserve(n);
  std::copy_if(dist.begin(), dist.end(), std::back_inserter(filtered_dist), [](float d) { return d >= 0.f; });
  dist = std::move(filtered_dist);

  // Here we start QMS computation
  qsm_start:

  // Now we can compute iter and clust.
  // What are iter and cluster varies depending on arbor versions
  // This is the part inspired from aRchi. We then cluster the point cloud
  // by iter and cluster
  std::vector<int>    iter = cut(dist, params.qsm.skeleton_node_distance);
  std::vector<int>   clust = cluster(wood, iter, params.qsm.dbscan_eps_distance);

  std::vector<std::pair<int, int>> iter_cluster;
  iter_cluster.reserve(n);
  for (std::size_t i = 0; i < n; ++i) iter_cluster.emplace_back(iter[i], clust[i]);

  // Build the QSM nodes and edges.
  // Returns the ids of the edges corresponding to each point so we can
  // later match back points to edges.
  auto ids = build_skeleton(wood, iter_cluster, params.qsm.max_d);

  // Assign edges IDs to the point cloud. Using the treeID attribute
  // because this is a free slot. Now each point has an ID corrsponding
  // to edges in the graph.
  for (size_t i = 0 ; i < wood.size() ; i ++) wood.set_treeid(i, ids[i]);
  ids.clear(); ids.shrink_to_fit();

  // Extremely rare case with so few points that we have no cluster
  // (seen once with a very bad DTM in Murray's data)
  if (graph.edges().size() == 0) {
    shift(tx, ty, tz);
    return;
  }

  // Is is still useful? build_skeleton() do it already?
  // TODO: review
  compute_topology();

  // Fix root issues (rare and likely even impossible with the new skeleton code)
  int n_root = count_nodes_connected_to_root();
  if (n_root == 0) throw std::runtime_error("Internal error in QSMbuilder::build. 0 root for this QSM. Please report.");
  if (n_root > 1) {
    ServiceLocator::logger()("Multiple nodes connected to root detected");
    fix_multiple_root();
  }

  // Compute branch order, axis IDs for the first time. We have 0 information on the edges
  // size. This is entirely based on the length of the branches in the graph
  compute_architecture(false);

  smooth_skeleton(params.qsm.smooth_steps);

  // Detect and fix some weird pattern in the bottom of the skeleton
  detect_weird_butt();

  // Measure the trees including
  // - initial measurements + quality estimation
  // - polynomial fitting
  // - reconstruction of missing radii
  bool success = construct_radii(wood, params.qsm.apex_radius);

  // A complete failure: we have 0 valid measurement, 0 radii. Retry one time with
  // less conservative parameter. This is very very rare.
  if (!success)
  {
    if (retried)
    {
      std::string msg = "[No valid measure] Not a single valid measure for this tree. The QSM is a pure reconstruction based on allometry";
      ServiceLocator::logger()("\033[33m" + msg + "\033[0m");
      graph.messages.push_back(msg);
    }
    else
    {
      std::string msg = "[No valid measure] Not a single valid measure for this tree. Retry with larger internode size";
      ServiceLocator::logger()("\033[33m" + msg + "\033[0m");
      retried = true;
      params.qsm.skeleton_node_distance *= 2;
      graph.clear();
      goto qsm_start;
    }
  }

  // Compute branch order, axis IDs for the second time. We have the radii so
  // now the computation is based on volume
  compute_architecture(true);

  smooth_skeleton(params.qsm.smooth_steps/2);

  // Special step if we confirm that the tree is likely broken
  // then the previous steps generated a very bad QSMs. Apply a special
  // reconstruction for dead trees.
  if (likely_broken && params.qsm.broken_detection_enabled)
  {
    refine_radii_broken(wood);
    compute_architecture(true);
    smooth_skeleton(params.qsm.smooth_steps/2);
  }

  // Prolongate the tree down to hag = 0
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

std::pair<int, int> which_min_max_z(const PointCloud& tree)
{
  double min_z = std::numeric_limits<double>::max();
  double max_z = -std::numeric_limits<double>::max();

  int min_index = -1;
  int max_index = -1;

  const std::size_t n = tree.size();

  for (std::size_t i = 0; i < n; ++i)
  {
    const double current_z = tree.get_z(i);

    if (current_z < min_z)
    {
      min_z = current_z;
      min_index = static_cast<int>(i);
    }

    if (current_z > max_z)
    {
      max_z = current_z;
      max_index = static_cast<int>(i);
    }
  }

  return {min_index, max_index};
}

PointCloud make_ground_seed(const PointCloud& wood)
{
  size_t n = wood.size();
  PointCloud gnd;
  double th = 0.05;
  while (gnd.size() == 0)
  {
    std::vector<bool> gnd_mask(n);
    for (size_t i = 0; i < n; i++)
    {
      if (wood.get_z(i) < th) //
        gnd_mask[i] = true;
    }
    gnd = wood.subset(gnd_mask);
    th += 0.05;
  }

  if (gnd.size() == wood.size() && wood.size() > 1)
    throw std::runtime_error("Internal error in QSMbuilder::build: gnd size == wood size. Please report.");

  return gnd;
}

PointCloud wood_subset(const PointCloud& tree, double z_threshold, int max_z_idx)
{
  std::vector<bool> wood_mask(tree.size(), false);
  for (std::size_t i = 0; i < tree.size(); ++i)
  {
    bool b = tree.is_wood(i) && (tree.get_z(i) >= z_threshold);
    wood_mask[i] = b;
  }
  wood_mask[max_z_idx] = true;

  return tree.subset(wood_mask);
}

}


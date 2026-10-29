/**
 * @file segment_semantic_api.cpp
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

#include "arbor.h"
#include "myomp.h"
#include "GraphBuilder.h"
#include "Grid3D.h"
#include "progress.h"

#include <limits>
#include <numeric>

namespace arbor::segment {

std::unique_ptr<Graph> build_semantic_graph(const PointCloud& core, const PointCloud& targets, const PointCloud& dtm, const settings::GraphParameters& params)
{
  Task task("Constructing the graph", {5, 2, 1, 2});

  GraphBuilder builder(params);

  task.step("Core layer");
  builder.add_core_layer(core);
  task.step("Target layer");
  builder.add_target_layer(core, targets);
  task.step("Seed layer");
  builder.add_seed_layer(core, dtm);
  task.step("Directed reachability");
  builder.fix_directed_reachability(core);
  return builder.get_graph();
}

std::vector<int> accumulate_passages(const PointCloud& core, const PointCloud& dtm, const settings::GraphParameters& params)
{
  if (core.size() == 0) throw std::runtime_error("accumulate_passages(): point cloud is empty.");
  if (dtm.size() == 0)  throw std::runtime_error("accumulate_passages(): seeds point cloud is empty.");

  Task task("Accumulating passages", {5, 5, 30, 25, 35});

  task.step("Decimating the point cloud");

  // Decimation
  std::vector<bool> keep1 = arbor::utils::homogeneization(core, params.decimation, true);
  PointCloud dec = core.subset(keep1, true);

  task.step("Discretizing scene space");

  std::vector<bool> keep2 = arbor::utils::homogeneization(core, params.space_res, false);
  PointCloud targets = core.subset(keep2, true);

  task.step("Constructing the graph");

  // Build graph
  std::unique_ptr<Graph> graph = build_semantic_graph(dec, targets, dtm, params);

  if (graph == nullptr) throw std::runtime_error("accumulate_passages: Failed to build graph (null pointer returned).");

  logging::debug("Graph size: " + Graph::format_bytes(graph->mem()));
  logging::debug("Graph nodes: " + std::to_string(graph->size()));

  size_t num_raw_points = core.size();
  size_t num_points = dec.size();

  // Targets are not graph nodes: each one hangs off its nearest core node with a zero cost
  Graph::NodeIDs target_anchors = std::move(graph->target_anchors);

  task.step("Graph resolution");

  // Precompute Dijkstra distances from the seeds, keeping the settle order
  Graph::NodeIDs settle_order;
  Graph::GraphCache cache = graph->compute_distances(&settle_order);

  // The graph is no longer needed, free it before the
  // passage accumulation to reduce peak memory.
  graph.reset();

  logging::debug("Graph cache size: " + Graph::format_bytes(Graph::cache_mem(cache)));

  task.step("Accumulating passages");

  Counter counter("Paths", settle_order.size());

  // The passage of a node is the number of target paths going through it, i.e. the number
  // of targets in its subtree of the shortest-path tree. Each target counts 1 on its anchor,
  // then counts are pushed to predecessors in reverse settle order, so a node's subtree is
  // complete before it is added to its parent. O(N) instead of one path walk per target.
  // Dropped targets and targets whose anchor is unreachable contribute nothing.
  const Graph::DistanceVector& distances = cache.first;
  const Graph::PredecessorVector& predecessors = cache.second;

  std::vector<int> passage(num_points, 0);
  for (Graph::NodeId anchor : target_anchors)
  {
    if (anchor != Graph::NO_NODE && distances[anchor] != std::numeric_limits<Graph::Cost>::infinity())
      passage[anchor] += 1;
  }

  std::size_t settled = 0;
  for (auto it = settle_order.rbegin(); it != settle_order.rend(); ++it)
  {
    counter.tick_at(settled++);
    if (counter.cancelled()) break;

    Graph::NodeId pred = predecessors[*it];
    if (pred != Graph::NO_NODE) passage[pred] += passage[*it];
  }

  counter.throw_if_cancelled();

  // Release the cache: both vectors need clear + shrink_to_fit to actually
  // return heap memory.
  cache.first.clear();
  cache.first.shrink_to_fit();
  cache.second.clear();
  cache.second.shrink_to_fit();
  Graph::NodeIDs().swap(settle_order);
  Graph::NodeIDs().swap(target_anchors);

  // Transfer passage values from dec back to core.
  // Points not in dec keep the default value of 0.
  std::vector<int> core_passage(num_raw_points, 0);

  size_t dec_idx = 0;
  for (size_t i = 0; i < num_raw_points; ++i)
  {
    if (keep1[i])
    {
      core_passage[i] = passage[dec_idx];
      ++dec_idx;
    }
  }

  return core_passage;
}

std::vector<bool> assign_wood_from_passage(const PointCloud& pc, const settings::SemanticParameters& params)
{
  if (pc.size() == 0)     throw std::runtime_error("assign_wood_from_passage(): point cloud is empty.");
  if (!pc.has_passage())  throw std::runtime_error("assign_wood_from_passage(): point cloud is missing required 'passage' attribute.");

  Task task("Pathfinder-based wood segmentation", {1, 2, 5});
  task.step("Selecting passage points");

  // Build a boolean mask of passage points and verify at least one exists.
  std::vector<bool> skeleton_mask(pc.size(), false);
  size_t passage_count = 0;
  for (size_t i = 0; i < pc.size(); ++i)
  {
    if (pc.get_passage(i) > params.min_passage)
    {
      skeleton_mask[i] = true;
      ++passage_count;
    }
  }

  if (passage_count == 0) throw std::runtime_error("assign_wood_from_passage: no passage points found (no points with passage > min_passage)");

  task.step("Building KD-tree");

  pc.build_index();

  task.step("k-nn search");

  Counter counter("k-nn search", pc.size());
  std::vector<bool> is_wood(pc.size(), false);

  const double dist_threshold_sq = params.wood_assignation_dist * params.wood_assignation_dist;
  const int k = params.wood_assignation_k;

  // Iterate directly over the skeleton mask, no need to copy passage points
  // into a separate PointCloud.
  #pragma omp parallel
  {
    std::vector<unsigned int> idx(k);
    std::vector<double>  dist(k);
    double q[3];

    #pragma omp for schedule(static)
    for (size_t i = 0; i < pc.size(); ++i)
    {
      if (counter.cancelled()) continue;
      counter.tick_at(i);

      if (!skeleton_mask[i]) continue;

      pc.get_point(i, q);
      pc.knn(q, k, idx, dist);

      for (int j = 0; j < k; ++j)
      {
        if (dist[j] < dist_threshold_sq)
          is_wood[idx[j]] = true;
      }
    }
  }

  counter.throw_if_cancelled();

  return is_wood;
}

std::vector<bool> assign_wood_from_high_likelihood(const PointCloud& pc, const settings::SemanticParameters& params)
{
  if (pc.size() == 0)     throw std::runtime_error("assign_wood_from_high_likelihood(): point cloud is empty.");
  if (!pc.has_foliage())  throw std::runtime_error("assign_wood_from_high_likelihood(): point cloud is missing required 'foliage' attribute.");
  if (!pc.has_pwood())    throw std::runtime_error("assign_wood_from_high_likelihood(): point cloud is missing required 'pwood' attribute.");

  Task task("High likelihood based wood segmentation", {1, 5, 1});
  task.step("Selecting high likelihood points");

  // Extract only high likelihood + already wood in previous step (assign_wood_from_passage)
  std::vector<bool> mask(pc.size(), false);
  for (size_t i = 0; i < pc.size(); ++i) {
    if (pc.is_wood(i) || pc.get_pwood(i) > params.high_pwood_threshold)
      mask[i] = true;
  }
  PointCloud wood = pc.subset(mask, true);

  task.step("Connected components");

  // Connected components to detect big clusters
  Grid3D grid(wood, params.connected_components_res);
  std::vector<int> cluster_ids = grid.connected_components(26);

  task.step("Connected components filtering");

  // Remove small clusters
  int max_id = *std::max_element(cluster_ids.begin(), cluster_ids.end());
  std::vector<int> counts(max_id + 1, 0);
  for (int id : cluster_ids) {
    if (id != 0)
      ++counts[id];
  }
  for (int& id : cluster_ids) {
    if (id != 0 && counts[id] < params.connected_components_min)
      id = 0;
  }

  // Map surviving cluster points back to the original point cloud
  std::vector<bool> is_wood(pc.size(), false);
  size_t j = 0;
  for (size_t i = 0; i < pc.size(); ++i)
  {
    if (mask[i]) {
      int cid = cluster_ids[j];
      if (cid > 0 && counts[cid] >= params.connected_components_min) {
        is_wood[i] = true;
      }

      ++j;
    }
  }

  return is_wood;
}

std::vector<bool> assign_wood_from_medium_likelihood(const PointCloud& pc, const settings::SemanticParameters& params)
{
  if (pc.size() == 0)     throw std::runtime_error("assign_wood_from_medium_likelihood(): point cloud is empty.");
  if (!pc.has_foliage())  throw std::runtime_error("assign_wood_from_medium_likelihood(): point cloud is missing required 'foliage' attribute.");
  if (!pc.has_pwood())    throw std::runtime_error("assign_wood_from_medium_likelihood(): point cloud is missing required 'pwood' attribute.");

  Task task("Medium likelihood based wood segmentation", {1, 4, 5, 1});
  task.step("Selecting medium likelihood points");

  // Extract only medium likelihood + already wood in previous steps
  std::vector<bool> mask(pc.size(), false);
  for (size_t i = 0; i < pc.size(); ++i) {
    if (pc.is_wood(i) || (pc.get_pwood(i) > params.medium_pwood_threshold && pc.get_pwood(i) < params.high_pwood_threshold))
      mask[i] = true;
  }
  PointCloud wood = pc.subset(mask, true);

  task.step("SOR noise filtering");

  // SOR. Detect noise
  std::vector<bool> is_noise = arbor::utils::sor(wood, params.medium_pwood_sor_k, params.medium_pwood_sor_m);

  // Remove noise
  is_noise.flip();
  wood = wood.subset(is_noise, true);
  is_noise.flip();

  task.step("Connected components");

  // Connected components to detect big clusters
  Grid3D grid(wood, params.connected_components_res);
  std::vector<int> cluster_ids = grid.connected_components(26);

  task.step("Connected components filtering");

  // Remove small clusters
  int max_id = *std::max_element(cluster_ids.begin(), cluster_ids.end());
  std::vector<int> counts(max_id + 1, 0);
  for (int id : cluster_ids) {
    if (id != 0)
      ++counts[id];
  }
  for (int& id : cluster_ids) {
    if (id != 0 && counts[id] < params.connected_components_min)
      id = 0;
  }

  // Map surviving cluster points back to the original point cloud.
  // j: index into the first subset (mask); k: index into the SOR-filtered subset.
  std::vector<bool> is_wood(pc.size(), false);
  size_t j = 0;
  size_t k = 0;

  for (size_t i = 0; i < pc.size(); ++i)
  {
    if (mask[i])
    {
      if (!is_noise[j])
      {
        int cid = cluster_ids[k];
        if (cid > 0 && counts[cid] >= params.connected_components_min)
          is_wood[i] = true;
        k++;
      }
      j++;
    }
  }

  return is_wood;
}

std::vector<bool> assign_wood_from_wood_dilatation(const PointCloud& pc, const settings::SemanticParameters& params)
{
  if (pc.size() == 0)     throw std::runtime_error("assign_wood_from_wood_dilatation: point cloud is empty.");
  if (!pc.has_foliage())  throw std::runtime_error("assign_wood_from_wood_dilatation: point cloud is missing required 'foliage' attribute.");

  // Dilate the current wood mask: neighbours of wood points within a radius
  // are also assigned as wood.  This recovers false negatives on thin branches.

  Task task("Dilatation based wood segmentation", {1, 2, 5});
  task.step("Selecting wood points");

  // Extract current wood points as a sub-cloud to iterate query coordinates
  std::vector<bool> is_wood(pc.size(), false);
  for (size_t i = 0; i < pc.size(); ++i)
  {
    if (pc.is_wood(i))
      is_wood[i] = true;
  }
  PointCloud wood = pc.subset(is_wood, true);

  task.step("Building KD-tree");

  pc.build_index();

  task.step("k-nn search");

  Counter counter("k-nn search", wood.size());

  const double max_dist_sq = params.wood_extra_reasignation_dist * params.wood_extra_reasignation_dist;
  const int k = params.wood_extra_reasignation_k;

  #pragma omp parallel
  {
    std::vector<unsigned int> idx(k);
    std::vector<double>  dist(k);
    double q[3];

    #pragma omp for schedule(static)
    for (size_t i = 0; i < wood.size(); ++i)
    {
      if (counter.cancelled()) continue;
      counter.tick_at(i);

      wood.get_point(i, q);
      pc.knn(q, k, idx, dist);

      for (int j = 0; j < k; ++j)
      {
        if (dist[j] <= max_dist_sq)
          is_wood[idx[j]] = true;
      }
    }
  }

  counter.throw_if_cancelled();

  return is_wood;
}

void segment_semantic(PointCloud& scene, const PointCloud& dtm, const settings::ArborParameters& par)
{
  if (scene.size() == 0)     throw std::runtime_error("segment_semantic: point cloud is empty.");
  if (!scene.has_foliage())  throw std::runtime_error("segment_semantic: point cloud is missing required 'foliage' attribute.");
  if (!scene.has_pwood())    throw std::runtime_error("segment_semantic: point cloud is missing required 'pwood' attribute.");
  if (!scene.has_userdata()) throw std::runtime_error("segment_semantic: point cloud is missing required 'UserData' attribute.");

  Task task("Semantic segmentation", {1, 60, 8, 8, 10, 8, 1});

  task.step("Partitioning");

  // The outputs start from scratch: a second run in place (e.g. with another cut height) gives the same
  // result as a first run
  for (size_t i = 0 ; i < scene.size() ; i++)
  {
    scene.set_userdata(i, static_cast<int>(PointCloud::UserData::TREE));
    scene.set_passage(i, 0);
    scene.set_foliage(i, 1);
  }

  // The low points are hidden at the end of the memory layout and will never be processed. The scene
  // gets its size and its order back when the partition goes out of scope, also on error or cancellation.
  ScopedPartition<PointCloud> partition(scene, [&](size_t i) { return scene.get_hag(i) > par.global.cut_above_ground; });
  size_t n = scene.size();
  if (n == 0) throw std::runtime_error("segment_semantic: no point above 'cut_above_ground'");

  // Flag low points with UserData = 1. Low points are beyond the size of the partitioned point cloud
  for (size_t i = n ; i < partition.total() ; i++) scene.set_userdata(i, static_cast<int>(PointCloud::UserData::LOW));

  task.step("Accumulating passages");
  std::vector<int> passages = accumulate_passages(scene, dtm, par.pathfinder);
  for (size_t i = 0 ; i < n ; i++) scene.set_passage(i, passages[i]);

  task.step("Pathfinder-based wood segmentation");
  std::vector<bool> path_finder_based_wood = assign_wood_from_passage(scene, par.semantic);
  for (size_t i = 0 ; i < n ; i++) scene.set_foliage(i, (int)!path_finder_based_wood[i]);

  task.step("High likelihood based wood segmentation");
  std::vector<bool> high_likelihood_based_wood   = assign_wood_from_high_likelihood(scene, par.semantic);
  task.step("Medium likelihood based wood segmentation");
  std::vector<bool> medium_likelihood_based_wood = assign_wood_from_medium_likelihood(scene, par.semantic);

  for (size_t i = 0 ; i < n ; i++)
  {
    bool wood = path_finder_based_wood[i] || high_likelihood_based_wood[i] || medium_likelihood_based_wood[i];
    scene.set_foliage(i, (int)(!wood));
  }

  task.step("Dilatation based wood segmentation");
  std::vector<bool> is_wood = assign_wood_from_wood_dilatation(scene, par.semantic);
  for (size_t i = 0 ; i < n ; i++) scene.set_foliage(i, (int)!is_wood[i]);

  task.step("Extra class 2 foliage re-assignation");
  for (size_t i = 0 ; i < n ; i++)
  {
    if (!scene.is_wood(i) && scene.get_pwood(i) > par.semantic.high_pwood_threshold)
      scene.set_foliage(i, 2);
  }
}

}

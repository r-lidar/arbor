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
#include "nanoflann.h"
#include "GraphBuilder.h"
#include "Grid3D.h"

#include <numeric>

using KDTree  = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 3>;
using index_t = nanoflann::KNNResultSet<double>::IndexType;

namespace arbor::segment {

Graph* build_semantic_graph(const PointCloud& core, const PointCloud& targets, const PointCloud& dtm, const settings::GraphParameters& params)
{
  GraphBuilder builder(params);

  builder.add_core_layer(core);
  builder.add_target_layer(core, targets);
  builder.add_seed_layer(core, dtm);
  builder.add_master_seed_layer();

  builder.fix_directed_reachability(core);

  return builder.get_graph();
}

std::vector<int> accumulate_passages(const PointCloud& core, const PointCloud& dtm, const settings::GraphParameters& params)
{
  if (core.size() == 0) throw std::runtime_error("accumulate_passages(): point cloud is empty.");
  if (dtm.size() == 0)  throw std::runtime_error("accumulate_passages(): seeds point cloud is empty.");

  ServiceLocator::logger()("Decimating the point cloud (1/10)");

  // Decimation
  std::vector<bool> keep1 = arbor::utils::homogeneization(core, params.decimation, true);
  PointCloud dec = core.subset(keep1, true);

  ServiceLocator::logger()("Discretizing scene space (2/10)");

  std::vector<bool> keep2 = arbor::utils::homogeneization(core, params.space_res, false);
  PointCloud targets = core.subset(keep2, true);

  ServiceLocator::logger()("Constructing the graph (3/10)");

  // Build graph
  Graph* graph = build_semantic_graph(dec, targets, dtm, params);

  ServiceLocator::logger()(" Graph size: " + Graph::format_bytes(graph->mem()));
  ServiceLocator::logger()(" Graph nodes: " + std::to_string(graph->adjacency_list.size()));

  if (graph == nullptr) throw std::runtime_error("segment_instance: Failed to build graph (null pointer returned).");

  size_t num_raw_points = core.size();
  size_t num_points = dec.size();
  size_t num_target = targets.size();
  size_t num_gnd    = dtm.size();

  std::vector<int> target_ids(num_target);
  std::vector<int> ground_ids(num_gnd);
  int master_id = num_points + num_target + num_gnd;

  std::iota(target_ids.begin(), target_ids.end(), num_points);
  std::iota(ground_ids.begin(), ground_ids.end(), num_points + num_target);

  // Global count vector
  std::vector<int> passage(num_points, 0);

  ServiceLocator::logger()("Graph resolution (4/10)");

  // Precompute distances for fast access
  Graph::GraphCache cache = graph->compute_distances(master_id);

  delete graph;

  ServiceLocator::logger()("  Graph cache size: " + Graph::format_bytes(Graph::cache_mem(cache)));

  ServiceLocator::logger()("Accumulating passages (5/10)");

  // Parallel loop over goal nodes
  #pragma omp parallel
  {
    std::vector<int> local_passage(num_points, 0);  // thread-local counts

    #pragma omp for schedule(dynamic, 100)
    for (size_t i = 0; i < num_target; ++i)
    {
      Graph::NodeId goal  = target_ids[i];

      auto [path, cost] = graph->findPath(master_id, goal, cache);

      for (size_t j = 0; j < path.size(); ++j)
      {
        Graph::NodeId id = path[j];
        if (id >= 0 && id < (int) num_points)
          local_passage[id] += 1;
      }
    }

    // Merge results into global passage safely
    #pragma omp critical
    {
      for (size_t i = 0; i < num_points; ++i)
        passage[i] += local_passage[i];
    }
  }

  // Release graph cache memory as it's no longer needed
  cache.first.clear();
  cache.first.shrink_to_fit();
  cache.second.clear();

  // Transfer passage values from dec back to core
  // Points not in dec get value 0
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

  ServiceLocator::logger()("Pathfinder-based wood segmentation (6/10)");

  // Filter pseudo-skeleton: points with passage > min_passage
  std::vector<bool> skeleton_mask(pc.size(), false);
  for (size_t i = 0; i < pc.size(); ++i)
  {
    if (pc.get_passage(i) > params.min_passage)
      skeleton_mask[i] = true;
  }
  PointCloud passages = pc.subset(skeleton_mask, true);

  if (passages.size() == 0) throw std::runtime_error("assign_wood_from_passage: no passage points found (no points with passage > min_passage)");

  ServiceLocator::logger()("  Building KDtree");

  // Spatial index of the skeleton
  KDTree tree(3, pc, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  tree.buildIndex();
  nanoflann::SearchParameters nanoparams;
  nanoparams.sorted = false;

  ServiceLocator::logger()("  k-nn search");

  // Each point close enough from a passage is assigned wood
  std::vector<bool> is_wood(pc.size(), false);

  // Precompute squared distance threshold to avoid repeated multiplication or sqrt
  const double dist_threshold_sq = params.wood_assignation_dist * params.wood_assignation_dist;
  const int k = params.wood_assignation_k;

  #pragma omp parallel
  {
    std::vector<index_t> idx(k);
    std::vector<double> dist(k);
    double q[3];

    #pragma omp for schedule(static)
    for (size_t i = 0; i < passages.size(); ++i)
    {
      passages.get_point(i, q);
      nanoflann::KNNResultSet<double> resultSet(k);
      resultSet.init(idx.data(), dist.data());
      tree.findNeighbors(resultSet, q, nanoparams);

      for (int j = 0 ; j < 10 ; j++)
      {
        if (dist[j] < dist_threshold_sq)
          is_wood[idx[j]] = true;
      }
    }
  }

  return is_wood;
}

std::vector<bool> assign_wood_from_high_likelihood(const PointCloud& pc, const settings::SemanticParameters& params)
{
  if (pc.size() == 0)     throw std::runtime_error("assign_wood_from_high_likelihood(): point cloud is empty.");
  if (!pc.has_foliage())  throw std::runtime_error("assign_wood_from_high_likelihood(): point cloud is missing required 'foliage' attribute.");
  if (!pc.has_pwood())    throw std::runtime_error("assign_wood_from_high_likelihood(): point cloud is missing required 'pwood' attribute.");

  ServiceLocator::logger()("High likelihood based wood segmentation (7/10)");

  // Extract only high likelihood + already wood in previous step (assign_wood_from_passage)
  std::vector<bool> mask(pc.size(), false);
  for (size_t i = 0; i < pc.size(); ++i) {
    if (pc.is_wood(i) || pc.get_pwood(i) > params.high_pwood_threshold)
      mask[i] = true;
  }
  PointCloud wood = pc.subset(mask, true);

  ServiceLocator::logger()("  Connected components computing");

  // Connected components to detect big clusters
  Grid3D grid(wood, params.connected_components_res);
  std::vector<int> cluster_ids = grid.connected_components(26);

  ServiceLocator::logger()("  Connected components filtering");

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

  // Assign original point cloud with wood/foliage
  std::vector<bool> is_wood(pc.size(), false);
  size_t j = 0;
  for (size_t i = 0; i < pc.size(); ++i)
  {
    if (mask[i]) {
      // In R, a point only survives if its cluster_id > 0
      // AND its count >= min.
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

  ServiceLocator::logger()("Medium likelihood based wood segmentation (8/10)");

  // Extract only medium likelihood + already wood in previous steps
  std::vector<bool> mask(pc.size(), false);
  for (size_t i = 0; i < pc.size(); ++i) {
    if (pc.is_wood(i) || (pc.get_pwood(i) > params.medium_pwood_threshold && pc.get_pwood(i) < params.high_pwood_threshold))
      mask[i] = true;
  }
  PointCloud wood = pc.subset(mask, true);

  ServiceLocator::logger()("  sor noise segmentation");

  // SOR. Detect noise
  std::vector<bool> is_noise = arbor::utils::sor(wood, params.medium_pwood_sor_k, params.medium_pwood_sor_m);

  // Remove noise
  is_noise.flip();
  wood = wood.subset(is_noise, true);
  is_noise.flip();

  ServiceLocator::logger()("  Connected components computing");

  // Connected components to detect big clusters
  Grid3D grid(wood, params.connected_components_res);
  std::vector<int> cluster_ids = grid.connected_components(26);

  ServiceLocator::logger()("  Connected components filtering");

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

  // Assign original point cloud with wood/foliage
  std::vector<bool> is_wood(pc.size(), false);
  size_t j = 0; // Index for wood_subset (mask)
  size_t k = 0; // Index for clustered_subset (survival_mask)

  for (size_t i = 0; i < pc.size(); ++i)
  {
    if (mask[i])
    {
      // Point was in the first subset. Did it survive SOR?
      if (!is_noise[j])
      {
        // Point was in the second subset. Did it survive Cluster Filtering?
        int cid = cluster_ids[k];
        if (cid > 0 && counts[cid] >= params.connected_components_min) {
          is_wood[i] = true;
        }
        k++; // Increment k only if the point survived SOR
      }
      j++; // Increment j every time mask[i] is true
    }
  }

  return is_wood;
}

std::vector<bool> assign_wood_from_wood_dilatation(const PointCloud& pc, const settings::SemanticParameters& params)
{
  if (pc.size() == 0)     throw std::runtime_error("assign_wood_from_medium_likelihood: point cloud is empty.");
  if (!pc.has_foliage())  throw std::runtime_error("assign_wood_from_medium_likelihood: point cloud is missing required 'foliage' attribute.");

  // We look at the neighboring points of the wood.  Points close to the wood
  // are wood points too. This assigns extra wood point is the branches and remove
  // some false negatives

  ServiceLocator::logger()("Dilatation based wood segmentation... (9/10)");

  // Extract wood points
  std::vector<bool> is_wood(pc.size(), false);
  for (size_t i = 0; i < pc.size(); ++i) {
    if (pc.is_wood(i))
      is_wood[i] = true;
  }
  PointCloud wood = pc.subset(is_wood, true);


  ServiceLocator::logger()("  Building KDtree");

  // Build KDTree on wood points
  KDTree tree(3, pc, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  tree.buildIndex();
  nanoflann::SearchParameters nanoparams;
  nanoparams.sorted = false;


  ServiceLocator::logger()("  knn search");

  const double max_dist_sq = params.wood_extra_reasignation_dist * params.wood_extra_reasignation_dist;
  const int k = params.wood_extra_reasignation_k;

  #pragma omp parallel
  {
    std::vector<index_t> idx(k);
    std::vector<double> dist(k);
    double q[3];

    #pragma omp for schedule(static)
    for (size_t i = 0; i < wood.size(); ++i)
    {
      wood.get_point(i, q);
      nanoflann::KNNResultSet<double> resultSet(k);
      resultSet.init(idx.data(), dist.data());
      tree.findNeighbors(resultSet, q, nanoparams);

      // If nearest wood point is within threshold distance, assign as wood
      for (int j = 0 ; j < 10 ; j++)
        if (dist[j] <= max_dist_sq)
          is_wood[idx[j]] = true;
    }
  }

  return is_wood;
}

void segment_semantic(PointCloud& scene, const PointCloud& dtm, const settings::ArborParameters& par)
{
  if (scene.size() == 0)     throw std::runtime_error("segment_semantic: point cloud is empty.");
  if (!scene.has_foliage())  throw std::runtime_error("segment_semantic: point cloud is missing required 'foliage' attribute.");
  if (!scene.has_pwood())    throw std::runtime_error("segment_semantic: point cloud is missing required 'pwood' attribute.");

  ServiceLocator::logger()("Partitioning...");
  scene.partition([&](size_t i) { return scene.get_hag(i) > par.global.cut_above_ground; });

  size_t n = scene.size();
  if (scene.size() == 0)     throw std::runtime_error("segment_semantic: no point above 'cut_above_ground'");

  std::vector<int> passages = accumulate_passages(scene, dtm, par.pathfinder);
  for (size_t i = 0 ; i < n ; i++) scene.set_passage(i, passages[i]);

  std::vector<bool> path_finder_based_wood = assign_wood_from_passage(scene, par.semantic);
  for (size_t i = 0 ; i < n ; i++) scene.set_foliage(i, (int)!path_finder_based_wood[i]);

  std::vector<bool> high_likelihood_based_wood = assign_wood_from_high_likelihood(scene, par.semantic);
  std::vector<bool> medium_likelihood_based_wood = assign_wood_from_medium_likelihood(scene, par.semantic);

  for (size_t i = 0 ; i < n ; i++) {
    bool wood = path_finder_based_wood[i] || high_likelihood_based_wood[i] || medium_likelihood_based_wood[i];
    scene.set_foliage(i, (int)(!wood));
  }

  std::vector<bool> is_wood = assign_wood_from_wood_dilatation(scene, par.semantic);
  for (size_t i = 0 ; i < n ; i++) scene.set_foliage(i, (int)!is_wood[i]);

  ServiceLocator::logger()("Extra class 2 foliage re-assignation... (10/10)");
  for (size_t i = 0 ; i < n ; i++)
  {
    if (!scene.is_wood(i) && scene.get_pwood(i) > par.semantic.high_pwood_threshold)
      scene.set_foliage(i, 2);
  }
}

}


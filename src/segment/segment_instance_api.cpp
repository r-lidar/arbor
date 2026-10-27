/**
 * @file segment_instance_api.cpp
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
#include "Grid3D.h"
#include "GraphBuilder.h"
#include "segment_overfitting.h"
#include "progress.h"

#include <numeric>
#include <queue>

using KDTree  = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 3>;
using index_t = nanoflann::KNNResultSet<double>::IndexType;

namespace arbor::segment {

std::unique_ptr<Graph> build_instance_graph(const PointCloud& core, const PointCloud& seeds, const settings::GraphParameters& params)
{
  if (core.size() == 0)     throw std::runtime_error("build_instance_graph: core point cloud is empty.");
  if (seeds.size() == 0)    throw std::runtime_error("build_instance_graph: seeds point cloud is empty.");
  if (!core.has_foliage())  throw std::runtime_error("build_instance_graph: core point cloud is missing required 'foliage' attribute.");
  if (!seeds.has_treeid())  throw std::runtime_error("build_instance_graph: seed point cloud is missing required 'treeid' attribute.");

  std::vector<bool> wood; wood.reserve(core.size());
  for (size_t i = 0; i < core.size(); ++i) wood.push_back(core.is_wood(i));

  Task task("Constructing the graph", {5, 1, 2});

  GraphBuilder builder(params);
  builder.set_wood(wood);
  task.step("Core layer");
  builder.add_core_layer(core);
  task.step("Seed layer");
  builder.add_seed_layer(core, seeds);
  task.step("Directed reachability");
  builder.fix_directed_reachability(core);
  return builder.get_graph();
}

void fix_small_isolated_low_clusters(PointCloud& las, double res = 0.05, int min_size = 200)
{
  const size_t n = las.size();

  // group indices by tree ID
  std::unordered_map<int, std::vector<int>> tree_to_indices;
  for (size_t i = 0; i < n; ++i)
  {
    if (!las.is_wood(i))        continue;
    if (las.get_hag(i) >= 3.0)  continue;
    int id = las.get_treeid(i);
    if (id <= 0)                continue; // NA from R or -1
    tree_to_indices[id].push_back(static_cast<int>(i));
  }

  Counter counter("Trees", tree_to_indices.size());

  // PER-TREE
  for (auto& [id, indices] : tree_to_indices)
  {
    counter.tick();
    counter.throw_if_cancelled();

    // Build a sub-cloud with Z scaled down by 0.1 to flatten the search
    PointCloud sub(indices.size());
    for (size_t k = 0; k < indices.size(); ++k)
    {
      int i = indices[k];
      sub.set_x(k, las.get_x(i));
      sub.set_y(k, las.get_y(i));
      sub.set_z(k, las.get_z(i) * 0.1);
    }

    // Run connected components on the flattened sub-cloud
    Grid3D grid(sub, res);
    std::vector<int> cluster_ids = grid.connected_components(26);

    // Count points per cluster
    std::unordered_map<int, int> counts;
    for (int cid : cluster_ids) counts[cid]++;

    if (counts.size() <= 1) continue;

    // Find the largest cluster
    int best_cluster = -1, best_count = -1;
    for (auto& [cid, cnt] : counts)
    {
      if (cnt > best_count) { best_count = cnt; best_cluster = cid; }
    }

    // Reclassify non-largest cluster points as foliage
    for (size_t k = 0; k < indices.size(); ++k)
    {
      if (cluster_ids[k] != best_cluster)
        las.set_foliage(indices[k], 1);
    }
  }
}

void segment_instance(PointCloud& core, const PointCloud& seeds, const settings::ArborParameters& params)
{
  if (core.size() == 0)     throw std::runtime_error("segment_instance: point cloud is empty.");
  if (seeds.size() == 0)    throw std::runtime_error("segment_instance: seeds point cloud contains 0 seed.");
  if (!core.has_foliage())  throw std::runtime_error("segment_instance: point cloud is missing required 'foliage' attribute.");
  if (!seeds.has_treeid())  throw std::runtime_error("segment_instance: seed point cloud is missing required 'treeid' attribute.");

  std::vector<double> weights = {1, 5, 30, 30, 10, 5};
  if (params.instance.oversegmentation_solver_enabled) weights.push_back(10);
  Task task("Instance segmentation", weights);

  task.step("Partitioning");
  core.partition([&](size_t i) { return core.get_hag(i) > params.global.cut_above_ground; });

  task.step("Decimating the point cloud");

  // Decimation
  std::vector<bool> keep = arbor::utils::homogeneization(core, params.pathfinder.decimation, true);
  PointCloud dec = core.subset(keep);

  size_t num_raw_pts = core.size();
  size_t num_points  = dec.size();

  task.step("Constructing the graph");

  // Build graph
  std::unique_ptr<Graph> graph = build_instance_graph(dec, seeds, params.pathfinder);

  if (graph == nullptr) throw std::runtime_error("segment_instance: Failed to build graph (null pointer returned).");

  logging::debug("Graph size: " + Graph::format_bytes(graph->mem()));

  task.step("Pathfinder");

  // Retrieve the closest seed for each point
  std::vector<double> distances;
  Graph::NodeIDs closest_seed;
  graph->shortest_paths_from_sources(distances, closest_seed);

  graph.reset();

  if (closest_seed.size() != num_points) throw std::runtime_error("segment_instance: Pathfinding returned incomplete results.");

  // Remap seed index -> seed id for each decimated point
  Graph::NodeIDs treeID(num_points, -1);
  for (size_t i = 0 ; i < num_points ; i++)
  {
    Graph::NodeId id = closest_seed[i];
    if (id != Graph::NO_NODE) treeID[i] = static_cast<Graph::NodeId>(seeds.get_treeid(id));
  }

  task.step("Assigning tree IDs to the dense point cloud");

  // Expand seed id from dec point cloud to core point cloud
  KDTree tree(3, dec, nanoflann::KDTreeSingleIndexAdaptorParams(10, nanoflann::KDTreeSingleIndexAdaptorFlags::SkipInitialBuildIndex));
  tree.buildIndex();
  nanoflann::SearchParameters nanoparams;
  nanoparams.sorted = false;

  Graph::NodeIDs ans(num_raw_pts);

  #pragma omp parallel
  {
    std::vector<index_t> idx(1);
    std::vector<double> dist(1);
    double q[3];

    #pragma omp for schedule(static)
    for (size_t i = 0; i < num_raw_pts; ++i)
    {
      core.get_point(i, q);
      nanoflann::KNNResultSet<double> resultSet(1);
      resultSet.init(idx.data(), dist.data());
      tree.findNeighbors(resultSet, q, nanoparams);
      index_t id = idx[0];
      ans[i] = treeID[id];
    }
  }

  for (size_t i = 0 ; i < core.size() ; i++)  core.set_treeid(i, ans[i]);

  task.step("Fixing low isolated wood clusters");
  fix_small_isolated_low_clusters(core);

  if (params.instance.oversegmentation_solver_enabled)
  {
    task.step("Overfitting enhancement");
    OverSegmentationResolver osr;
    osr.run(core);
  }
}


std::vector<float> dist2root(const PointCloud& core, const PointCloud& dtm, const settings::GraphParameters& params)
{
  if (core.size() == 0) throw std::runtime_error("dist2root(): point cloud is empty.");
  if (dtm.size() == 0)  throw std::runtime_error("dist2root(): seeds point cloud is empty.");

  size_t num_points = core.size();

  Task task("Distance to root", {10, 10, 1});

  task.step("Constructing the graph");

  GraphBuilder builder(params);
  builder.add_core_layer(core);
  builder.add_seed_layer(core, dtm);
  builder.fix_directed_reachability2(core);

  std::unique_ptr<Graph> graph = builder.get_graph();
  if (graph == nullptr)
    throw std::runtime_error("dist2root(): failed to build graph (null pointer).");

  // Run Dijkstra from the seeds
  task.step("Pathfinder");

  Graph::NodeIDs settle_order;
  Graph::GraphCache cache = graph->compute_distances(&settle_order);
  const auto& graph_distances = cache.first;
  const auto& predecessors = cache.second;

  graph.reset();

  logging::debug("Graph cache size: " + Graph::format_bytes(Graph::cache_mem(cache)));

  // For each core point, walk the predecessor chain and accumulate Euclidean distance
  task.step("Computing euclidean distance to ground");

  std::vector<float> euclidean_distance_to_root(num_points, -1.0);

  const Graph::NodeId NO_PRED = Graph::NO_NODE;

  // Precompute per-node euclidean contribution (from node -> its predecessor) O(N) point fetches total
  std::vector<float> edge_euclidean(num_points, 0.0);
  for (size_t i = 0; i < num_points; ++i)
  {
    if (graph_distances[i] == std::numeric_limits<Graph::Cost>::infinity()) continue;

    // Nodes without predecessor are connected to a seed: they are roots
    Graph::NodeId prev = predecessors[i];
    if (prev == NO_PRED) continue;

    double c[3], p[3];
    core.get_point(i, c);
    core.get_point(prev, p);
    const double dx = c[0]-p[0], dy = c[1]-p[1], dz = c[2]-p[2];
    edge_euclidean[i] = static_cast<float>(std::sqrt(dx*dx + dy*dy + dz*dz));
  }

  // Accumulate along the tree using memoization
  // Process in Dijkstra settle order: a topological order of the shortest-path tree, so a
  // node's predecessor is always processed before the node itself. Only reached nodes are in it.
  for (Graph::NodeId idx : settle_order)
  {
    Graph::NodeId prev = predecessors[idx];

    if (prev == NO_PRED)
    {
      euclidean_distance_to_root[idx] = 0.0; // root itself
      continue;
    }

    euclidean_distance_to_root[idx] = euclidean_distance_to_root[prev] + edge_euclidean[idx];
  }

  return euclidean_distance_to_root;
}

}

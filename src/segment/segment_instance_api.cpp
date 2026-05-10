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

#include <chrono>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <queue>

using KDTree  = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 3>;
using index_t = nanoflann::KNNResultSet<double>::IndexType;

namespace arbor::segment {

Graph* build_instance_graph(const PointCloud& core, const PointCloud& seeds, const settings::GraphParameters& params)
{
  if (core.size() == 0)     throw std::runtime_error("build_instance_graph: core point cloud is empty.");
  if (seeds.size() == 0)    throw std::runtime_error("build_instance_graph: seeds point cloud is empty.");
  if (!core.has_foliage())  throw std::runtime_error("build_instance_graph: core point cloud is missing required 'foliage' attribute.");
  if (!seeds.has_treeid())  throw std::runtime_error("build_instance_graph: seed point cloud is missing required 'treeid' attribute.");

  std::vector<bool> wood; wood.reserve(core.size());
  for (size_t i = 0; i < core.size(); ++i) wood.push_back(core.is_wood(i));

  GraphBuilder builder(params);
  builder.set_wood(wood);

  builder.add_core_layer(core);
  builder.add_seed_layer(core, seeds);
  builder.add_master_seed_layer();

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

  // PER-TREE
  for (auto& [id, indices] : tree_to_indices)
  {
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

  ServiceLocator::logger()("Partitioning...");
  core.partition([&](size_t i) { return core.get_hag(i) > params.global.cut_above_ground; });

  const auto t0 = std::chrono::steady_clock::now();

  ServiceLocator::logger()("Instance segmentation start");
  ServiceLocator::logger()("Decimating the point cloud... (1/4)");

  // Decimation
  std::vector<bool> keep = arbor::utils::homogeneization(core, params.pathfinder.decimation, true);
  PointCloud dec = core.subset(keep);

  size_t num_raw_pts = core.size();
  size_t num_points  = dec.size();
  size_t num_seeds   = seeds.size();

  ServiceLocator::logger()("Constructing the graph (2/4)");

  // Build graph
  Graph* graph = build_instance_graph(dec, seeds, params.pathfinder);

  if (graph == nullptr) throw std::runtime_error("segment_instance: Failed to build graph (null pointer returned).");

  // Indexes of the seeds in the graph
  Graph::NodeIDs seeds_ids(num_seeds);
  std::iota(seeds_ids.begin(), seeds_ids.end(), num_points);

  ServiceLocator::logger()("Pathfinder (3/4)");

  // Retrieve the closest seed for each point
  std::vector<double> distances;
  Graph::NodeIDs closest_nodeids;
  graph->shortest_paths_from_node(seeds_ids, distances, closest_nodeids);
  delete graph;

  if (closest_nodeids.size() != num_points+num_seeds+1) throw std::runtime_error("segment_instance: Pathfinding returned incomplete results.");

  // Remap seed index -> seed id for each decimated point
  Graph::NodeIDs treeID(num_points, -1);
  for (size_t i = 0 ; i < num_points ; i++)
  {
    Graph::NodeId id = closest_nodeids[i];
    if (id != -1) treeID[i] = static_cast<Graph::NodeId>(seeds.get_treeid(id-num_points));
  }

  ServiceLocator::logger()("Assigning tree IDs to the dense point cloud (4/4)");

  // Expand seed id from dec point cloud to core point cloud
  KDTree tree(3, dec, nanoflann::KDTreeSingleIndexAdaptorParams(10));
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

  const auto t1 = std::chrono::steady_clock::now();
  const std::chrono::duration<double> elapsed = t1 - t0;

  std::ostringstream oss;
  oss << std::fixed << std::setprecision(1) << elapsed.count();

  ServiceLocator::logger()("Fix low isolated wood clusters");
  fix_small_isolated_low_clusters(core);

  ServiceLocator::logger()("Instance segmentation completed in " + oss.str() + " s");
}


std::vector<float> dist2root(const PointCloud& core, const PointCloud& dtm, const settings::GraphParameters& params)
{
  ServiceLocator::logger()("Computing shortest paths to ground");

  if (core.size() == 0) throw std::runtime_error("dist2root(): point cloud is empty.");
  if (dtm.size() == 0)  throw std::runtime_error("dist2root(): seeds point cloud is empty.");

  ServiceLocator::logger()(" Constructing the graph");
  GraphBuilder builder(params);
  builder.add_core_layer(core);
  builder.add_seed_layer(core, dtm);
  builder.add_master_seed_layer();

  Graph* graph = builder.get_graph();
  if (graph == nullptr) throw std::runtime_error("dist2root(): failed to build graph (null pointer).");

  const size_t num_points = core.size();
  const size_t num_gnd    = dtm.size();
  const Graph::NodeId master_id = static_cast<Graph::NodeId>(num_points + num_gnd);

  // Force-connect the highest core point to its nearest neighbors
  // The graph is directed and distance-limited, so the apex can end up isolated.
  // We guarantee it has edges regardless of builder constraints.
  {
    constexpr size_t K = 10;

    // Find the apex (max Z).
    size_t apex_id = 0;
    double apex_z  = -std::numeric_limits<double>::max();
    for (size_t i = 0; i < num_points; ++i)
    {
      double p[3];
      core.get_point(i, p);
      if (p[2] > apex_z) { apex_z = p[2]; apex_id = i; }
    }

    // Collect the K nearest neighbours with a max-heap of size K.
    using Entry = std::pair<double, size_t>;
    std::priority_queue<Entry> heap;
    double ac[3];
    core.get_point(apex_id, ac);

    for (size_t i = 0; i < num_points; ++i)
    {
      if (i == apex_id) continue;
      double p[3];
      core.get_point(i, p);
      const double dx = ac[0]-p[0], dy = ac[1]-p[1], dz = ac[2]-p[2];
      double d  = std::pow(std::sqrt(dx*dx + dy*dy + dz*dz), params.power);
      d = std::pow(d, params.power);

      heap.push({d, i});
      if (heap.size() > K) heap.pop(); // evict the farthest
    }

    // Add directional edges for each neighbor.
    while (!heap.empty())
    {
      const auto [d, nb_id] = heap.top(); heap.pop();
      graph->add_edge(static_cast<Graph::NodeId>(nb_id), static_cast<Graph::NodeId>(apex_id), d);
      //graph->add_edge(static_cast<Graph::NodeId>(apex_id), static_cast<Graph::NodeId>(nb_id), d);
    }
  }


  // Run Dijkstra from master seed
  ServiceLocator::logger()(" Dijkstra");
  Graph::GraphCache cache = graph->compute_distances(master_id);
  const auto& graph_distances = cache.first;
  const auto& predecessors = cache.second;

  // Identify reachable and isolated sets
  ServiceLocator::logger()(" Identifying isolated sub-graph");
  std::vector<Graph::NodeId> reachable, isolated;

  for (size_t i = 0; i < num_points; ++i)
  {
    if (graph_distances[i] == std::numeric_limits<Graph::Cost>::infinity())
      isolated.push_back(static_cast<Graph::NodeId>(i));
    else
      reachable.push_back(static_cast<Graph::NodeId>(i));
  }

  // Add fallback edges: isolated -> nearest reachable core point
  if (isolated.size() > 100 && !reachable.empty())
  {
    ServiceLocator::logger()("\033[33m " + std::to_string(isolated.size()) + " isolated points detected: restart and expand.\033[0m");
    for (size_t ii = 0; ii < isolated.size(); ii += 2)
    {
      Graph::NodeId iso = isolated[ii];
      double ci[3];
      core.get_point(iso, ci);

      double best_dist = std::numeric_limits<double>::max();
      Graph::NodeId best_node = reachable[0];

      for (size_t k = 0; k < reachable.size(); k += 2)
      {
        Graph::NodeId reach = reachable[k];
        double cr[3];
        core.get_point(reach, cr);
        double dx = ci[0]-cr[0], dy = ci[1]-cr[1], dz = ci[2]-cr[2];
        double d  = std::sqrt(dx*dx + dy*dy + dz*dz);
        d = std::pow(d, params.power);
        if (d < best_dist) { best_dist = d; best_node = reach; }
      }

      graph->add_edge(best_node, iso, best_dist);
    }

    // Second Dijkstra run on augmented graph
    cache = graph->compute_distances(master_id);
  }

  delete graph;

  // For each core point, walk the predecessor chain and accumulate Euclidean distance
  ServiceLocator::logger()(" Computing euclidian distance to ground");
  std::vector<float> euclidean_distance_to_root(num_points, -1.0);

  const Graph::NodeId NO_PRED = -1;

  // Precompute per-node euclidean contribution (from node -> its predecessor)
  // O(N) point fetches total
  std::vector<float> edge_euclidean(num_points, 0.0);
  for (size_t i = 0; i < num_points; ++i)
  {
    if (graph_distances[i] == std::numeric_limits<Graph::Cost>::infinity()) continue;
    if (i >= predecessors.size()) continue;
    Graph::NodeId prev = predecessors[i];
    if (prev == NO_PRED) continue;
    if (prev < static_cast<Graph::NodeId>(num_points))
    {
      double c[3], p[3];
      core.get_point(i, c);
      core.get_point(prev, p);
      const double dx = c[0]-p[0], dy = c[1]-p[1], dz = c[2]-p[2];
      edge_euclidean[i] = static_cast<float>(std::sqrt(dx*dx + dy*dy + dz*dz));
    }
  }

  // Accumulate along the tree using memoization
  // Process in topological order (by graph_distances = BFS/Dijkstra order)
  std::vector<size_t> order(num_points);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b)
  {
    return graph_distances[a] < graph_distances[b];
  });
  for (size_t idx : order)
  {
    if (graph_distances[idx] == std::numeric_limits<Graph::Cost>::infinity()) continue;
    Graph::NodeId prev = (idx < predecessors.size()) ? predecessors[idx] : NO_PRED;
    if (prev == NO_PRED)
    {
      euclidean_distance_to_root[idx] = 0.0; // root itself
      continue;
    }
    double parent_dist = (prev < static_cast<Graph::NodeId>(num_points)) ? euclidean_distance_to_root[prev] : 0.0;
    euclidean_distance_to_root[idx] = parent_dist + edge_euclidean[idx];
  }

  return euclidean_distance_to_root;
}

}

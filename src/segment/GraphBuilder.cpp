/**
 * @file GraphBuilder.cpp
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

#include "GraphBuilder.h"
#include "PointCloud.h"
#include "myomp.h"

#include <vector>
#include <unordered_map>
#include <cmath>

namespace arbor::segment {

GraphBuilder::GraphBuilder(const settings::GraphParameters& p)
{
  params = p;
  set_angle_penalty(p.angle_penalty);  // validate
  graph = new Graph();
}

GraphBuilder::~GraphBuilder() { if (graph_owner) delete graph; }
Graph* GraphBuilder::get_graph() { graph_owner = false; return graph; }
void GraphBuilder::set_wood(const std::vector<bool>& x) { wood = x; }

void GraphBuilder::set_angle_penalty(const std::vector<float>& x)
{
  constexpr std::size_t expected_size = 181;

  if (x.size() != expected_size)
  {
    throw std::runtime_error(
        "Invalid angle penalty factor vector size: expected " +
          std::to_string(expected_size) +
          ", got " +
          std::to_string(x.size()) + "."
    );
  }

  params.angle_penalty = x;
}


// ---------------------------------------------------------
// 1. Core Layer (semi bidirectional)
// ---------------------------------------------------------

void GraphBuilder::add_core_layer(const PointCloud& core)
{
  if (total_core_nodes > 0)   throw std::runtime_error("Core layer already populated");
  if (total_target_nodes > 0) throw std::runtime_error("Core layer must be populated first");
  if (total_seed_nodes > 0)   throw std::runtime_error("Core layer must be populated first");
  if (total_master_nodes > 0) throw std::runtime_error("Core layer must be populated first");

  // +1 because nanoflann's kNN includes the query point itself
  params.k++;

  int n_points = core.size();
  bool use_wood = wood.size() > 0;

  offset_points    = 0;
  total_core_nodes = n_points;
  total_nodes      = n_points;

  graph->ensure_size(total_nodes);
  graph->reserve_edges(params.k);

  // Build the KD-tree once; it is reused by add_target_layer, add_seed_layer,
  // and fix_directed_reachability, then released at the end of that last call.
  auto start = std::chrono::high_resolution_clock::now();
  core_index = std::make_unique<KDTree>(3, core, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  core_index->buildIndex();
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  ServiceLocator::logger()(" KDTree creation : " + std::to_string(duration.count()) + " s");

  #pragma omp parallel
  {
    std::vector<size_t> idx(params.k);
    std::vector<double> dist(params.k);

    // Thread-local storage for edges
    std::vector<std::tuple<Graph::NodeId, Graph::NodeId, Graph::Cost>> local_edges;

    // For each point we connect to its knn. The current point is 'from'.
    #pragma omp for schedule(static)
    for (int from = 0; from < n_points; ++from)
    {
      double q[3];
      core.get_point(from, q);

      nanoflann::KNNResultSet<double> result(params.k);
      result.init(&idx[0], &dist[0]);
      core_index->findNeighbors(result, q, nanoflann::SearchParameters());

      // For each neighbour, compute the cost to connect from -> to
      for (int j = 0; j < params.k; ++j)
      {
        int to = idx[j];
        if (to == from) continue;               // skip self (0-nn)
        float cost = std::sqrt(dist[j]);
        if (cost > params.max_gap) continue;    // gap too large: no connection
        cost = std::pow(cost, params.power);

        // q already holds coord_from ,  no second fetch needed
        double coord_to[3];
        core.get_point(to, coord_to);

        // Direction penalty: moving upward is cheap, downward is expensive.
        float dx = q[0] - coord_to[0];
        float dy = q[1] - coord_to[1];
        float dz = q[2] - coord_to[2];
        float magnitude = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (magnitude < 1e-12f) continue;
        float cos_theta = -dz / magnitude;
        if (params.downward) cos_theta = -cos_theta;
        float angle_deg = std::acos(std::clamp(cos_theta, -1.0f, 1.0f)) * 180.0f / M_PI;
        int   angle     = std::round(angle_deg);
        cost *= params.angle_penalty[angle];  // use the rounded integer index

        // Optional wood/foliage cost modifiers
        if (use_wood)
        {
          bool is_wood1 = wood[from];
          bool is_wood2 = wood[to];
          if (is_wood1 && is_wood2)   cost *= params.wood2wood;
          else if (!is_wood1 && !is_wood2) cost *= params.leaf2leaf;
          else if (is_wood1 && !is_wood2)  cost *= params.wood2leaf;
        }

        local_edges.emplace_back(from, to, cost);
      }
    }

    // Parallel reduction: merge thread-local edges into the shared graph
    #pragma omp critical
    {
      for (auto& e : local_edges)
        graph->add_edge(std::get<0>(e), std::get<1>(e), std::get<2>(e));
    }
  }
}

// ---------------------------------------------------------
// 2. Target Layer (core point -> target node)
// ---------------------------------------------------------

// Each target point is connected to its nearest core point with cost 0.

void GraphBuilder::add_target_layer(const PointCloud& core, const PointCloud& target)
{
  if (total_target_nodes > 0) throw std::runtime_error("Target layer already populated");
  if (total_core_nodes == 0)  throw std::runtime_error("Target layer must be populated after core layer");
  if (total_seed_nodes > 0)   throw std::runtime_error("Target layer must be populated before seed layer");
  if (total_master_nodes > 0) throw std::runtime_error("Target layer must be populated before master layer");

  int n_target = target.size();

  offset_targets     = total_nodes;
  total_target_nodes = n_target;     // fixed: was incorrectly set to core.size()
  total_nodes       += n_target;

  graph->ensure_size(total_nodes);
  graph->reserve_edges(params.k);

  // Reuse the KD-tree built in add_core_layer
  std::vector<size_t> idx(params.k);
  std::vector<double> dist(params.k);

  for (int i = 0; i < n_target; ++i)
  {
    double q[3];
    target.get_point(i, q);
    nanoflann::KNNResultSet<double> result(params.k);
    result.init(&idx[0], &dist[0]);
    core_index->findNeighbors(result, q, nanoflann::SearchParameters());

    // Zero-cost edge from the nearest core point to this target node
    int from = idx[0];
    int to   = i + offset_targets;
    graph->add_edge(from, to, 0.0f);
  }
}

// ---------------------------------------------------------
// 3. Seed Layer (seed -> core points)
// ---------------------------------------------------------

void GraphBuilder::add_seed_layer(const PointCloud& core, const PointCloud& seeds)
{
  if (total_seed_nodes > 0)   throw std::runtime_error("Seed layer already populated");
  if (total_core_nodes == 0)  throw std::runtime_error("Seed layer must be populated after core layer");
  if (total_master_nodes > 0) throw std::runtime_error("Seed layer must be populated before master layer");

  int k        = params.k_seed;
  int n_points = seeds.size();

  offset_seeds     = total_nodes;
  total_seed_nodes = n_points;
  total_nodes     += n_points;

  graph->ensure_size(total_nodes);
  graph->reserve_edges(k);

  // Reuse the KD-tree built in add_core_layer
  std::vector<size_t> idx(k);
  std::vector<double> dist(k);

  for (int i = 0; i < n_points; ++i)
  {
    double q[3];
    seeds.get_point(i, q);
    nanoflann::KNNResultSet<double> result(k);
    result.init(&idx[0], &dist[0]);
    core_index->findNeighbors(result, q, nanoflann::SearchParameters());

    for (int j = 0; j < k; ++j)
    {
      float cost = std::sqrt(dist[j]);
      cost = std::pow(cost, params.power);
      int from = i + offset_seeds;
      int to   = idx[j];
      graph->add_edge(from, to, cost);
    }
  }
}

// ---------------------------------------------------------
// 4. Master Seed Layer (master -> all seeds, cost 0)
// ---------------------------------------------------------

void GraphBuilder::add_master_seed_layer()
{
  if (total_master_nodes > 0) throw std::runtime_error("Master layer already populated");
  if (total_core_nodes == 0)  throw std::runtime_error("Master layer must be populated after core layer");
  if (total_seed_nodes == 0)  throw std::runtime_error("Seed layer must be populated before master layer");

  offset_master      = total_nodes;
  total_master_nodes = 1;
  total_nodes       += 1;

  graph->ensure_size(total_nodes);
  graph->reserve_edges(total_seed_nodes);

  for (int i = 0; i < total_seed_nodes; ++i)
    graph->add_edge(offset_master, offset_seeds + i, 0.0f);
}

// ---------------------------------------------------------
// 5. Fix gap and unreachable points
// ---------------------------------------------------------

void GraphBuilder::fix_directed_reachability(const PointCloud& cloud)
{
  ServiceLocator::logger()(" Test directed reachability");

  const Graph::NodeId source = get_range_master().second;
  const size_t N = static_cast<size_t>(total_core_nodes);
  if (N < 2)
  {
    core_index.reset();
    return;
  }

  // DFS from master seed: returns core node IDs unreachable from source
  auto find_isolated = [&]() -> std::vector<Graph::NodeId>
  {
    const size_t total = graph->adjacency_list.size();
    std::vector<bool> reachable(total, false);
    std::vector<Graph::NodeId> stack;
    stack.reserve(total);
    stack.push_back(source);
    reachable[source] = true;

    while (!stack.empty())
    {
      Graph::NodeId u = stack.back(); stack.pop_back();
      for (const auto& e : graph->adjacency_list[u])
      {
        if (e.destination >= 0 && !reachable[e.destination])
        {
          reachable[e.destination] = true; stack.push_back(e.destination);
        }
      }
    }

    std::vector<Graph::NodeId> isolated;
    for (size_t i = 0; i < N; ++i)
    {
      if (!reachable[i])
        isolated.push_back(static_cast<Graph::NodeId>(i));
    }

    return isolated;
  };

  std::vector<Graph::NodeId> isolated_ids = find_isolated();
  if (isolated_ids.empty())
  {
    core_index.reset();  // no isolated nodes, free tree memory early
    return;
  }

  ServiceLocator::logger()(" " + std::to_string(isolated_ids.size()) + " isolated nodes: rebuilding bidirectionnal edges with relaxed params");

  // Relaxed parameters for the retry pass
  const int   k_retry   = params.k * 4;
  const float gap_retry = params.max_gap;// * 4.0f;
  const bool  use_wood  = !wood.empty();

  // Reuse the KD-tree built in add_core_layer (same cloud, same parameters)
  std::vector<std::tuple<Graph::NodeId, Graph::NodeId, Graph::Cost>> new_edges;
  new_edges.reserve(isolated_ids.size() * k_retry * 2);

  std::vector<size_t> idx(k_retry);
  std::vector<double> sq_dist(k_retry);

  for (Graph::NodeId from : isolated_ids)
  {
    // Wipe old edges, they were confined within the isolated sub-graph
    graph->adjacency_list[from].clear();

    double q[3];
    cloud.get_point(from, q);

    nanoflann::KNNResultSet<double> result(k_retry);
    result.init(idx.data(), sq_dist.data());
    core_index->findNeighbors(result, q, nanoflann::SearchParameters());

    for (int j = 0; j < k_retry; ++j)
    {
      int to = static_cast<int>(idx[j]);
      if (to == from) continue;

      float eucl = std::sqrt(sq_dist[j]);
      if (eucl > gap_retry) continue;

      // q already holds coord_from ,  no second fetch needed
      double coord_to[3];
      cloud.get_point(to, coord_to);

      float dx        = q[0] - coord_to[0];
      float dy        = q[1] - coord_to[1];
      float dz        = q[2] - coord_to[2];
      float magnitude = std::sqrt(dx*dx + dy*dy + dz*dz);
      if (magnitude < 1e-12f) continue;

      float base_cost = std::pow(eucl, params.power);

      // Forward edge: from -> to
      {
        float cos_theta = -dz / magnitude;
        if (params.downward) cos_theta = -cos_theta;
        int angle = std::round(std::acos(std::clamp(cos_theta, -1.0f, 1.0f)) * 180.0f / M_PI);
        float cost = base_cost * params.angle_penalty[angle];

        if (use_wood)
        {
          bool w1 = wood[from], w2 = wood[to];
          if      ( w1 &&  w2) cost *= params.wood2wood;
          else if (!w1 && !w2) cost *= params.leaf2leaf;
          else if ( w1 && !w2) cost *= params.wood2leaf;
        }

        new_edges.emplace_back(from, to, cost);
      }

      // Reverse edge: to -> from (direction flipped, so dz and wood roles swap)
      {
        float cos_theta = dz / magnitude;   // +dz because direction is reversed
        if (params.downward) cos_theta = -cos_theta;
        int angle = std::round(std::acos(std::clamp(cos_theta, -1.0f, 1.0f)) * 180.0f / M_PI);
        float cost = base_cost * params.angle_penalty[angle];

        if (use_wood)
        {
          bool w1 = wood[to], w2 = wood[from];   // roles swapped
          if      ( w1 &&  w2) cost *= params.wood2wood;
          else if (!w1 && !w2) cost *= params.leaf2leaf;
          else if ( w1 && !w2) cost *= params.wood2leaf;
        }

        new_edges.emplace_back(to, from, cost);
      }
    }
  }

  for (auto& [a, b, c] : new_edges)
    graph->add_edge(a, b, c);

  // Verify
  std::vector<Graph::NodeId> still_isolated = find_isolated();
  if (still_isolated.empty())
    ServiceLocator::logger()(" All isolated nodes resolved");
  else
    ServiceLocator::logger()(" " + std::to_string(still_isolated.size()) + " nodes still isolated after retry (genuine gaps)");

  // KD-tree is no longer needed ,  release its memory before returning
  core_index.reset();
}

void GraphBuilder::fix_directed_reachability2(const PointCloud& cloud)
{
  ServiceLocator::logger()(" Test directed reachability");

  core_index.reset();

  const Graph::NodeId source = get_range_master().second;
  const size_t N = static_cast<size_t>(total_core_nodes);
  if (N < 2) return;

  // ---------------------------------------------------------------
  // Step 1: BFS from source following directed edges.
  //         Mark every node reachable from the master seed.
  // ---------------------------------------------------------------
  const size_t total = graph->adjacency_list.size();
  std::vector<bool> reachable(total, false);
  {
    std::vector<Graph::NodeId> stack;
    stack.reserve(total);
    stack.push_back(source);
    reachable[source] = true;
    while (!stack.empty())
    {
      Graph::NodeId u = stack.back(); stack.pop_back();
      for (const auto& e : graph->adjacency_list[u])
        if (e.destination >= 0 && !reachable[e.destination])
        { reachable[e.destination] = true; stack.push_back(e.destination); }
    }
  }

  // ---------------------------------------------------------------
  // Step 2: Partition core nodes into reachable vs. isolated.
  // ---------------------------------------------------------------
  std::vector<Graph::NodeId> isolated_ids;
  std::vector<Graph::NodeId> reachable_ids;
  isolated_ids.reserve(N);
  reachable_ids.reserve(N);

  for (size_t i = 0; i < N; ++i)
  {
    if (reachable[i])
      reachable_ids.push_back(static_cast<Graph::NodeId>(i));
    else
      isolated_ids.push_back(static_cast<Graph::NodeId>(i));
  }

  if (isolated_ids.empty()) return;

  // ---------------------------------------------------------------
  // Step 3: Find connected components within the isolated subgraph.
  //
  // We do BFS/DFS following only edges whose *destination* is also
  // isolated.  Because fix_directed_reachability() may have already added
  // bidirectional edges inside isolated sub-graphs, following
  // outgoing edges is sufficient to discover all cluster members.
  //
  // iso_component[node_id] = cluster index, or -1 when not isolated.
  // This gives O(1) neighbour lookup without a hash map.
  // ---------------------------------------------------------------
  std::vector<int> iso_component(N, -1);   // -1 = reachable (not isolated)
  for (Graph::NodeId id : isolated_ids)
    iso_component[id] = INT_MAX;           // isolated but unvisited sentinel

  int num_clusters = 0;

  for (Graph::NodeId seed : isolated_ids)
  {
    if (iso_component[seed] != INT_MAX) continue;  // already labelled

    // BFS over the isolated subgraph
    std::vector<Graph::NodeId> queue = { seed };
    iso_component[seed] = num_clusters;

    for (size_t qi = 0; qi < queue.size(); ++qi)
    {
      Graph::NodeId u = queue[qi];
      for (const auto& e : graph->adjacency_list[u])
      {
        Graph::NodeId v = e.destination;
        if (v >= 0 && v < static_cast<Graph::NodeId>(N) &&
            iso_component[v] == INT_MAX)
        {
          iso_component[v] = num_clusters;
          queue.push_back(v);
        }
      }
    }

    ++num_clusters;
  }

  ServiceLocator::logger()(" " +
    std::to_string(isolated_ids.size()) +
    " isolated nodes grouped into " +
    std::to_string(num_clusters) +
    " clusters -> bridging");

  // ---------------------------------------------------------------
  // Step 4: Build a KD-tree over the reachable core points only.
  // ---------------------------------------------------------------
  PointCloud reachable_cloud = cloud.subset(reachable_ids, true);
  KDTree tree(3, reachable_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  tree.buildIndex();
  nanoflann::SearchParameters search_params;
  search_params.sorted = true;

  // ---------------------------------------------------------------
  // Step 5: For every isolated node, find its nearest reachable
  //         neighbour and track the best (shortest) bridge per
  //         cluster.  One bridge edge is all that is needed: once
  //         the entry point is reachable, BFS propagates through
  //         the cluster's own internal edges.
  // ---------------------------------------------------------------
  struct Bridge
  {
    Graph::NodeId iso_node   = -1;
    Graph::NodeId reach_node = -1;
    double        sq_dist    = std::numeric_limits<double>::infinity();
  };
  std::vector<Bridge> best(num_clusters);

  std::vector<size_t> idx(1);
  std::vector<double> sq_dists(1);

  for (Graph::NodeId iso : isolated_ids)
  {
    double q[3];
    cloud.get_point(iso, q);

    nanoflann::KNNResultSet<double> result_set(1);
    result_set.init(idx.data(), sq_dists.data());
    tree.findNeighbors(result_set, q, search_params);

    int comp = iso_component[iso];
    if (sq_dists[0] < best[comp].sq_dist)
    {
      best[comp].sq_dist    = sq_dists[0];
      best[comp].iso_node   = iso;
      best[comp].reach_node = reachable_ids[idx[0]];
    }
  }

  // ---------------------------------------------------------------
  // Step 6: Commit one bidirectional bridge per cluster.
  //         No max_gap guard here: we unconditionally connect
  //         because the first pass already tried relaxed params
  //         and failed. The graph must be connected.
  // ---------------------------------------------------------------
  int bridges_added = 0;
  for (const Bridge& b : best)
  {
    if (b.iso_node == -1) continue;
    double d = std::sqrt(b.sq_dist);
    const Graph::Cost cost = static_cast<Graph::Cost>(std::pow(d, params.power));
    graph->add_edge(b.reach_node, b.iso_node, cost);
    graph->add_edge(b.iso_node,   b.reach_node, cost);
    ++bridges_added;
  }

  // ---------------------------------------------------------------
  // Step 7: Verify, one final reachability scan.
  // ---------------------------------------------------------------
  {
    std::fill(reachable.begin(), reachable.end(), false);
    std::vector<Graph::NodeId> stack;
    stack.push_back(source);
    reachable[source] = true;
    while (!stack.empty())
    {
      Graph::NodeId u = stack.back(); stack.pop_back();
      for (const auto& e : graph->adjacency_list[u])
      {
        if (e.destination >= 0 && !reachable[e.destination])
        {
          reachable[e.destination] = true; stack.push_back(e.destination);
        }
      }
    }

    size_t still_isolated = 0;
    for (size_t i = 0; i < N; ++i)
    {
      if (!reachable[i])
        ++still_isolated;
    }

    if (still_isolated == 0)
      ServiceLocator::logger()(" All isolated clusters resolved (" +  std::to_string(bridges_added) + " bridge edges added)");
    else
      ServiceLocator::logger()(" " + std::to_string(still_isolated) + " nodes still unreachable after pass 2 (genuine gaps)");
  }
}

int GraphBuilder::get_num_cores()   const { return total_core_nodes; }
int GraphBuilder::get_num_targets() const { return total_target_nodes; }
int GraphBuilder::get_num_seeds()   const { return total_seed_nodes; }
int GraphBuilder::get_num_master()  const { return total_master_nodes; }
std::pair<int, int> GraphBuilder::get_range_core()    const { return {offset_points,  offset_points  + total_core_nodes   - 1}; }
std::pair<int, int> GraphBuilder::get_range_targets() const { return {offset_targets, offset_targets + total_target_nodes - 1}; }
std::pair<int, int> GraphBuilder::get_range_seed()    const { return {offset_seeds,   offset_seeds   + total_seed_nodes   - 1}; }
std::pair<int, int> GraphBuilder::get_range_master()  const { return {offset_master,  offset_master  + total_master_nodes - 1}; }

}

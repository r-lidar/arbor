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
#include "progress.h"

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>

namespace arbor::segment {

GraphBuilder::GraphBuilder(const settings::GraphParameters& p)
{
  params = p;
  set_angle_penalty(p.angle_penalty);  // validate
}

GraphBuilder::~GraphBuilder() = default;
std::unique_ptr<Graph> GraphBuilder::get_graph() { return std::move(graph); }
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
  if (graph != nullptr) throw std::runtime_error("Core layer already populated");

  // +1 because nanoflann's kNN includes the query point itself
  params.k++;

  int n_points = core.size();
  bool use_wood = wood.size() > 0;

  total_core_nodes = n_points;

  // At most k-1 edges per node: the query point itself is skipped. It can only be missing from
  // its own kNN if k other points share its exact position, and those are skipped as well.
  const int stride = std::max(1, params.k - 1);
  graph = std::make_unique<Graph>(n_points, stride);

  core.build_index(); // no op if already built

  #pragma omp parallel
  {
    std::vector<unsigned int> idx(params.k);
    std::vector<double> dist(params.k);

    // For each point we connect to its knn. The current point is 'from'. Each thread writes
    // only the rows of the points it processes, so no synchronization is needed.
    #pragma omp for schedule(static)
    for (int from = 0; from < n_points; ++from)
    {
      double q[3];
      core.get_point(from, q);
      core.knn(q, params.k, idx, dist);

      Graph::Node* row = graph->row(from);
      int degree = 0;

      // For each neighbour, compute the cost to connect from -> to
      for (int j = 0; j < params.k && degree < stride; ++j)
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

        row[degree++] = {static_cast<Graph::NodeId>(to), cost};
      }

      // The slots are not initialized by Graph: mark the unused ones
      for (; degree < stride; ++degree)
        row[degree] = {Graph::NO_NODE, 0.0f};
    }
  }
}

// ---------------------------------------------------------
// 2. Target Layer (core point -> target)
// ---------------------------------------------------------

// Each target point hangs off its nearest core point with cost 0. Targets are leaves, so they
// are not stored as nodes: only the nearest core point (anchor) is recorded.

void GraphBuilder::add_target_layer(const PointCloud& core, const PointCloud& target)
{
  if (total_target_nodes > 0) throw std::runtime_error("Target layer already populated");
  if (graph == nullptr)       throw std::runtime_error("Target layer must be populated after core layer");

  int n_target = target.size();
  total_target_nodes = n_target;

  core.build_index(); // no op if already built

  Graph::NodeIDs& anchors = graph->target_anchors;
  anchors.assign(n_target, Graph::NO_NODE);

  #pragma omp parallel
  {
    std::vector<unsigned int> idx(1);
    std::vector<double> dist(1);

    #pragma omp for schedule(static)
    for (int i = 0; i < n_target; ++i)
    {
      double q[3];
      target.get_point(i, q);
      core.knn(q, 1, idx, dist);
      anchors[i] = idx[0];
    }
  }
}

// ---------------------------------------------------------
// 3. Seed Layer (seed -> core points)
// ---------------------------------------------------------

// Seeds are not stored as nodes either. With the master seed connected to every seed at cost 0,
// a shortest path starts at the seed offering the cheapest connection to its first core node.
// So each core node reached by a seed becomes a Dijkstra source with the cheapest seed cost.

void GraphBuilder::add_seed_layer(const PointCloud& core, const PointCloud& seeds)
{
  if (total_seed_nodes > 0) throw std::runtime_error("Seed layer already populated");
  if (graph == nullptr)     throw std::runtime_error("Seed layer must be populated after core layer");

  const int k        = params.k_seed;
  const int n_points = seeds.size();
  const size_t n_core = total_core_nodes;

  total_seed_nodes = n_points;

  core.build_index(); // no op if already built

  std::vector<Graph::Cost> best_cost(n_core, std::numeric_limits<Graph::Cost>::infinity());
  Graph::NodeIDs best_seed(n_core, Graph::NO_NODE);

  // kNN queries run in parallel by chunks, then the candidates are reduced sequentially in seed
  // order so that ties are resolved exactly as in a sequential search (first seed wins).
  const int chunk = 4096;
  std::vector<unsigned int> chunk_idx(static_cast<size_t>(chunk) * k);
  std::vector<Graph::Cost> chunk_cost(static_cast<size_t>(chunk) * k);

  for (int start = 0; start < n_points; start += chunk)
  {
    const int end = std::min(start + chunk, n_points);

    #pragma omp parallel
    {
      std::vector<unsigned int> idx(k);
      std::vector<double> dist(k);

      #pragma omp for schedule(static)
      for (int i = start; i < end; ++i)
      {
        double q[3];
        seeds.get_point(i, q);
        core.knn(q, k, idx, dist);

        size_t offset = static_cast<size_t>(i - start) * k;
        for (int j = 0; j < k; ++j)
        {
          float cost = std::sqrt(dist[j]);
          cost = std::pow(cost, params.power);
          chunk_idx[offset + j]  = idx[j];
          chunk_cost[offset + j] = cost;
        }
      }
    }

    for (int i = start; i < end; ++i)
    {
      size_t offset = static_cast<size_t>(i - start) * k;
      for (int j = 0; j < k; ++j)
      {
        unsigned int to = chunk_idx[offset + j];
        if (chunk_cost[offset + j] < best_cost[to])
        {
          best_cost[to] = chunk_cost[offset + j];
          best_seed[to] = i;
        }
      }
    }
  }

  graph->sources.clear();
  for (size_t i = 0; i < n_core; ++i)
  {
    if (best_seed[i] != Graph::NO_NODE)
      graph->sources.push_back({static_cast<Graph::NodeId>(i), best_cost[i], best_seed[i]});
  }
}

// ---------------------------------------------------------
// 4. Fix gap and unreachable points
// ---------------------------------------------------------

std::vector<bool> GraphBuilder::reachable_from_sources() const
{
  const size_t N = graph->size();
  std::vector<bool> reachable(N, false);
  std::vector<Graph::NodeId> stack;
  stack.reserve(N);

  for (const Graph::Source& s : graph->sources)
  {
    if (!reachable[s.node])
    {
      reachable[s.node] = true;
      stack.push_back(s.node);
    }
  }

  while (!stack.empty())
  {
    Graph::NodeId u = stack.back(); stack.pop_back();
    graph->for_each_edge(u, [&](const Graph::Node& e)
    {
      if (!reachable[e.destination])
      {
        reachable[e.destination] = true;
        stack.push_back(e.destination);
      }
    });
  }

  return reachable;
}

void GraphBuilder::fix_directed_reachability(const PointCloud& cloud)
{
  if (graph == nullptr)       throw std::runtime_error("Reachability must be fixed after core layer");
  if (total_seed_nodes == 0)  throw std::runtime_error("Reachability must be fixed after seed layer");

  const size_t N = static_cast<size_t>(total_core_nodes);

  if (N < 2) return;

  // Traversal from the sources: returns core node IDs unreachable from them
  auto find_isolated = [&]() -> std::vector<Graph::NodeId>
  {
    std::vector<bool> reachable = reachable_from_sources();

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
    return;
  }

  logging::debug(std::to_string(isolated_ids.size()) + " isolated nodes: rebuilding bidirectionnal edges with relaxed params");

  // Relaxed parameters for the retry pass
  const int   k_retry   = params.k * 4;
  const float gap_retry = params.max_gap;// * 4.0f;
  const bool  use_wood  = !wood.empty();

  // Reuse the KD-tree built in add_core_layer (same cloud, same parameters)
  std::vector<Graph::Edge> new_edges;
  new_edges.reserve(isolated_ids.size() * k_retry * 2);

  cloud.build_index();

  std::vector<unsigned int> idx(k_retry);
  std::vector<double> sq_dist(k_retry);

  std::vector<bool> wiped(N, false);

  for (Graph::NodeId from : isolated_ids)
  {
    // Wipe old edges, they were confined within the isolated sub-graph
    graph->clear_row(from);
    wiped[from] = true;

    double q[3];
    cloud.get_point(from, q);
    cloud.knn(q, k_retry, idx, sq_dist);

    for (int j = 0; j < k_retry; ++j)
    {
      unsigned int to = idx[j];
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

        new_edges.push_back({from, {to, cost}});
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

        new_edges.push_back({to, {from, cost}});
      }
    }
  }

  graph->add_extra_edges(std::move(new_edges));

  // When targets were stored as nodes, the zero-cost edge to a target was in the out-edges of
  // its core node and was wiped with them. Drop those targets to keep the same results.
  for (Graph::NodeId& anchor : graph->target_anchors)
  {
    if (anchor != Graph::NO_NODE && wiped[anchor])
      anchor = Graph::NO_NODE;
  }

  // Verify
  std::vector<Graph::NodeId> still_isolated = find_isolated();
  if (still_isolated.empty())
    logging::debug("All isolated nodes resolved");
  else
    logging::info(std::to_string(still_isolated.size()) + " nodes still isolated after retry (genuine gaps)");
}

void GraphBuilder::fix_directed_reachability2(const PointCloud& cloud)
{
  if (graph == nullptr)       throw std::runtime_error("Reachability must be fixed after core layer");
  if (total_seed_nodes == 0)  throw std::runtime_error("Reachability must be fixed after seed layer");

  const size_t N = static_cast<size_t>(total_core_nodes);
  if (N < 2) return;

  // ---------------------------------------------------------------
  // Step 1: Traversal from the sources following directed edges.
  //         Mark every node reachable from the seeds.
  // ---------------------------------------------------------------
  std::vector<bool> reachable = reachable_from_sources();

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
      graph->for_each_edge(u, [&](const Graph::Node& e)
      {
        Graph::NodeId v = e.destination;
        if (iso_component[v] == INT_MAX)
        {
          iso_component[v] = num_clusters;
          queue.push_back(v);
        }
      });
    }

    ++num_clusters;
  }

  logging::debug(
    std::to_string(isolated_ids.size()) +
    " isolated nodes grouped into " +
    std::to_string(num_clusters) +
    " clusters -> bridging");

  // ---------------------------------------------------------------
  // Step 4: Build a KD-tree over the reachable core points only.
  // ---------------------------------------------------------------
  PointCloud reachable_cloud = cloud.subset(reachable_ids, Attr::None);
  reachable_cloud.build_index();

  // ---------------------------------------------------------------
  // Step 5: For every isolated node, find its nearest reachable
  //         neighbour and track the best (shortest) bridge per
  //         cluster.  One bridge edge is all that is needed: once
  //         the entry point is reachable, BFS propagates through
  //         the cluster's own internal edges.
  // ---------------------------------------------------------------
  struct Bridge
  {
    Graph::NodeId iso_node   = Graph::NO_NODE;
    Graph::NodeId reach_node = Graph::NO_NODE;
    double        sq_dist    = std::numeric_limits<double>::infinity();
  };
  std::vector<Bridge> best(num_clusters);

  std::vector<unsigned int> idx(1);
  std::vector<double> sq_dists(1);

  for (Graph::NodeId iso : isolated_ids)
  {
    double q[3];
    cloud.get_point(iso, q);
    reachable_cloud.knn(q, 1, idx, sq_dists);

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
  std::vector<Graph::Edge> bridges;
  bridges.reserve(2 * best.size());
  for (const Bridge& b : best)
  {
    if (b.iso_node == Graph::NO_NODE) continue;
    double d = std::sqrt(b.sq_dist);
    const Graph::Cost cost = static_cast<Graph::Cost>(std::pow(d, params.power));
    bridges.push_back({b.reach_node, {b.iso_node, cost}});
    bridges.push_back({b.iso_node,   {b.reach_node, cost}});
  }
  int bridges_added = bridges.size() / 2;
  graph->add_extra_edges(std::move(bridges));

  // ---------------------------------------------------------------
  // Step 7: Verify, one final reachability scan.
  // ---------------------------------------------------------------
  {
    reachable = reachable_from_sources();

    size_t still_isolated = 0;
    for (size_t i = 0; i < N; ++i)
    {
      if (!reachable[i])
        ++still_isolated;
    }

    if (still_isolated == 0)
      logging::debug("All isolated clusters resolved (" +  std::to_string(bridges_added) + " bridge edges added)");
    else
      logging::info(std::to_string(still_isolated) + " nodes still unreachable after pass 2 (genuine gaps)");
  }
}

}

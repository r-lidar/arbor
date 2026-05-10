/**
 * @file Graph.cpp
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

#include <limits>
#include <algorithm>
#include <iostream>
#include <queue>

#include "Graph.h"

namespace arbor::segment {

void Graph::add_edge(NodeId source, NodeId destination, Cost cost)
{
  adjacency_list[source].push_back({destination, cost});
}

void Graph::ensure_size(size_t n)
{
  if (adjacency_list.size() < n)
    adjacency_list.resize(n);
}

void Graph::reserve_edges(NodeId node, size_t capacity)
{
  if (node >= 0 && node < static_cast<NodeId>(adjacency_list.size()))
    adjacency_list[node].reserve(capacity);
}

std::pair<Graph::DistanceVector, Graph::PredecessorVector> Graph::compute_distances(NodeId start) const
{
  DistanceVector distances(adjacency_list.size(), std::numeric_limits<Cost>::infinity());
  distances[start] = 0.0f;

  PredecessorVector predecessors(adjacency_list.size(), -1);

  using QueueNode = std::pair<Cost, NodeId>;
  std::priority_queue<QueueNode, std::vector<QueueNode>, std::greater<>> open_set;
  open_set.push({0.0f, start});

  while (!open_set.empty())
  {
    auto [current_cost, current] = open_set.top();
    open_set.pop();

    if (current_cost > distances[current]) continue;

    for (const Node& neighbor : adjacency_list[current])
    {
      Cost new_cost = current_cost + neighbor.cost;
      if (new_cost < distances[neighbor.destination])
      {
        distances[neighbor.destination] = new_cost;
        predecessors[neighbor.destination] = current;
        open_set.push({new_cost, neighbor.destination});
      }
    }
  }

  return {distances, predecessors};
}

// --- Path reconstruction ---
std::pair<Graph::Path, Graph::Cost> Graph::findPath(NodeId start, NodeId goal, const std::pair<DistanceVector, PredecessorVector>& precomputed_data) const
{
  const auto& [distances, predecessors] = precomputed_data;

  if (distances[goal] == std::numeric_limits<Cost>::infinity())
    return {{}, -1.0f};

  Path path;
  for (NodeId node = goal; node != start; )
  {
    if (node < 0 || node >= static_cast<NodeId>(predecessors.size()) || predecessors[node] == -1)
      return {{}, -1.0f};
    path.push_back(node);
    node = predecessors[node];
  }

  path.push_back(start);
  std::reverse(path.begin(), path.end());
  return {path, distances[goal]};
}

void Graph::shortest_paths_from_node(const NodeIDs& nodes, std::vector<double>& distances, NodeIDs& closest_nodeids) const
{
  const double INF = std::numeric_limits<double>::infinity();
  size_t N = adjacency_list.size();

  distances.assign(N, INF);
  closest_nodeids.assign(N, -1); // invalid default

  using PQElement = std::pair<double, NodeId>;
  std::priority_queue<PQElement, std::vector<PQElement>, std::greater<>> pq;

  // Initialize ground nodes
  for (NodeId g : nodes)
  {
    distances[g] = 0.0;
    closest_nodeids[g] = g;
    pq.push({0.0, g});
  }

  size_t processed = 0;
  size_t next_report = N / 20; // report every 5%
  if (next_report == 0) next_report = 1; // avoid division by zero

  while (!pq.empty())
  {
    auto [dist_u, u] = pq.top();
    pq.pop();

    if (dist_u > distances[u]) continue; // outdated entry

    // Progress tracking
    ++processed;
    if (processed % next_report == 0)
    {
      //double pct = 100.0 * processed / N;
      //std::cout << "Progress: " << static_cast<int>(pct) << "% (" << processed << "/" << N << " nodes processed)\r";
      //std::cout.flush();
    }

    for (const auto& e : adjacency_list[u])
    {
      NodeId v = e.destination;
      double new_dist = dist_u + e.cost;
      if (new_dist < distances[v])
      {
        distances[v] = new_dist;
        closest_nodeids[v] = closest_nodeids[u];
        pq.push({new_dist, v});
      }
    }
  }

  //std::cout << "\nDone. Processed " << processed << " nodes.\n";
}

}


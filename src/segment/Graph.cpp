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

// ---------------------------------------------------------------------------
// ensure_size / add_edge
// ---------------------------------------------------------------------------

void Graph::ensure_size(size_t n)
{
  if (n > num_nodes_)
    num_nodes_ = n;
}

void Graph::add_edge(NodeId source, NodeId destination, Cost cost)
{
  if (!csr_valid_)
  {
    // Build phase: accumulate in the flat COO buffer.
    coo_src_.push_back(source);
    coo_edges_.push_back({destination, cost});
  }
  else
  {
    // CSR already built: place in the small overflow map so Dijkstra can
    // still find these edges without a full rebuild.
    overflow_[source].push_back({destination, cost});
  }
}

// ---------------------------------------------------------------------------
// finalize — convert COO → CSR, then free the COO buffers.
// ---------------------------------------------------------------------------

void Graph::finalize() const
{
  if (csr_valid_) return;

  const size_t E = coo_src_.size();

  // 1. Count the out-degree of every node.
  row_ptr_.assign(num_nodes_ + 1, 0);
  for (NodeId src : coo_src_)
    row_ptr_[static_cast<size_t>(src) + 1]++;

  // 2. Exclusive prefix-sum → row_ptr_[i] is the start offset for node i.
  for (size_t i = 1; i <= num_nodes_; ++i)
    row_ptr_[i] += row_ptr_[i - 1];

  // 3. Fill adj_data_ using a temporary position array (a copy of row_ptr_).
  adj_data_.resize(E);
  {
    std::vector<size_t> pos(row_ptr_.begin(), row_ptr_.end());
    for (size_t i = 0; i < E; ++i)
    {
      const NodeId src = coo_src_[i];
      adj_data_[pos[static_cast<size_t>(src)]++] = coo_edges_[i];
    }
  } // pos freed here

  // 4. Release the (now redundant) COO buffers.
  coo_src_   = std::vector<NodeId>();
  coo_edges_ = std::vector<Node>();

  csr_valid_ = true;
}

// ---------------------------------------------------------------------------
// compute_distances — Dijkstra from a single source.
// ---------------------------------------------------------------------------

std::pair<Graph::DistanceVector, Graph::PredecessorMap> Graph::compute_distances(NodeId start) const
{
  finalize();

  DistanceVector distances(num_nodes_, std::numeric_limits<Cost>::infinity());
  distances[start] = 0.0f;

  PredecessorMap predecessors;

  using QueueNode = std::pair<Cost, NodeId>;
  std::priority_queue<QueueNode, std::vector<QueueNode>, std::greater<>> open_set;
  open_set.push({0.0f, start});

  while (!open_set.empty())
  {
    auto [current_cost, current] = open_set.top();
    open_set.pop();

    if (current_cost > distances[current]) continue;

    // CSR neighbours
    for (size_t i = row_ptr_[current]; i < row_ptr_[current + 1]; ++i)
    {
      const Node& neighbor = adj_data_[i];
      Cost new_cost = current_cost + neighbor.cost;
      if (new_cost < distances[neighbor.destination])
      {
        distances[neighbor.destination] = new_cost;
        predecessors[neighbor.destination] = current;
        open_set.push({new_cost, neighbor.destination});
      }
    }

    // Overflow neighbours (edges added after CSR was built — rare)
    auto ov = overflow_.find(current);
    if (ov != overflow_.end())
    {
      for (const Node& neighbor : ov->second)
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
  }

  return {distances, predecessors};
}

// ---------------------------------------------------------------------------
// findPath — path reconstruction from precomputed Dijkstra data.
// ---------------------------------------------------------------------------

std::pair<Graph::Path, Graph::Cost> Graph::findPath(NodeId start, NodeId goal, const std::pair<DistanceVector, PredecessorMap>& precomputed_data) const
{
  const auto& [distances, predecessors] = precomputed_data;

  if (distances[goal] == std::numeric_limits<Cost>::infinity())
    return {{}, -1.0f};

  Path path;
  for (NodeId node = goal; node != start; )
  {
    auto it = predecessors.find(node);
    if (it == predecessors.end())
      return {{}, -1.0f};
    path.push_back(node);
    node = it->second;
  }

  path.push_back(start);
  std::reverse(path.begin(), path.end());
  return {path, distances[goal]};
}

// ---------------------------------------------------------------------------
// shortest_paths_from_node — multi-source Dijkstra.
// ---------------------------------------------------------------------------

void Graph::shortest_paths_from_node(const NodeIDs& nodes, std::vector<double>& distances, NodeIDs& closest_nodeids) const
{
  finalize();

  const double INF = std::numeric_limits<double>::infinity();

  distances.assign(num_nodes_, INF);
  closest_nodeids.assign(num_nodes_, -1);

  using PQElement = std::pair<double, NodeId>;
  std::priority_queue<PQElement, std::vector<PQElement>, std::greater<>> pq;

  for (NodeId g : nodes)
  {
    distances[g] = 0.0;
    closest_nodeids[g] = g;
    pq.push({0.0, g});
  }

  size_t processed = 0;
  size_t next_report = num_nodes_ / 20;
  if (next_report == 0) next_report = 1;

  while (!pq.empty())
  {
    auto [dist_u, u] = pq.top();
    pq.pop();

    if (dist_u > distances[u]) continue;

    ++processed;
    if (processed % next_report == 0)
    {
      //double pct = 100.0 * processed / num_nodes_;
      //std::cout << "Progress: " << static_cast<int>(pct) << "% (" << processed << "/" << num_nodes_ << " nodes processed)\r";
      //std::cout.flush();
    }

    // CSR neighbours
    for (size_t i = row_ptr_[u]; i < row_ptr_[u + 1]; ++i)
    {
      const Node& e = adj_data_[i];
      double new_dist = dist_u + e.cost;
      if (new_dist < distances[e.destination])
      {
        distances[e.destination] = new_dist;
        closest_nodeids[e.destination] = closest_nodeids[u];
        pq.push({new_dist, e.destination});
      }
    }

    // Overflow neighbours (edges added after CSR was built — rare)
    auto ov = overflow_.find(u);
    if (ov != overflow_.end())
    {
      for (const Node& e : ov->second)
      {
        double new_dist = dist_u + e.cost;
        if (new_dist < distances[e.destination])
        {
          distances[e.destination] = new_dist;
          closest_nodeids[e.destination] = closest_nodeids[u];
          pq.push({new_dist, e.destination});
        }
      }
    }
  }

  //std::cout << "\nDone. Processed " << processed << " nodes.\n";
}

}

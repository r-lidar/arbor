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

#include <iomanip>
#include <limits>
#include <algorithm>
#include <queue>
#include <sstream>
#include <stdexcept>

#include "Graph.h"
#include "progress.h"

namespace arbor::segment {

Graph::Graph(size_t num_nodes, size_t stride) :
  num_nodes_(num_nodes),
  stride_(stride),
  slots_(new Node[num_nodes * stride]) // default-initialized: no serial pre-fill pass
{
}

void Graph::clear_row(NodeId u)
{
  Node* r = row(u);
  std::fill(r, r + stride_, Node{NO_NODE, 0.0f});
}

void Graph::add_extra_edges(std::vector<Edge> edges)
{
  if (edges.empty()) return;

  // Counting sort into CSR: no sort and no copy of the input. Existing edges of a node stay
  // first, new ones follow in insertion order.
  const size_t old_total = extra_edges_.size();
  const size_t total = old_total + edges.size();
  if (total > std::numeric_limits<uint32_t>::max())
    throw std::runtime_error("Graph: too many extra edges.");

  // offset[u + 1] = number of edges of u, then prefix sum: offset[u + 1] = end of the edges of u
  std::vector<uint32_t> offset(num_nodes_ + 1, 0);
  if (!extra_offset_.empty())
  {
    for (size_t u = 0; u < num_nodes_; ++u)
      offset[u + 1] = extra_offset_[u + 1] - extra_offset_[u];
  }
  for (const Edge& e : edges) offset[e.source + 1]++;
  for (size_t u = 0; u < num_nodes_; ++u) offset[u + 1] += offset[u];

  // Fill from the end of each node's range: new edges in reverse order, then existing ones.
  // offset[u + 1] is used as a cursor and ends up at the start of the range of u + 1.
  std::vector<Node> placed(total);
  for (size_t i = edges.size(); i-- > 0; )
    placed[--offset[edges[i].source + 1]] = edges[i].node;

  std::vector<Edge>().swap(edges);

  if (!extra_offset_.empty())
  {
    for (size_t u = 0; u < num_nodes_; ++u)
    {
      for (uint32_t e = extra_offset_[u + 1]; e-- > extra_offset_[u]; )
        placed[--offset[u + 1]] = extra_edges_[e];
    }
  }

  // The cursors have moved one node back: offset[u + 1] is now the start of u, shift them
  for (size_t u = 0; u < num_nodes_; ++u) offset[u] = offset[u + 1];
  offset[num_nodes_] = static_cast<uint32_t>(total);

  extra_offset_ = std::move(offset);
  extra_edges_ = std::move(placed);
}

Graph::GraphCache Graph::compute_distances(NodeIDs* settle_order) const
{
  DistanceVector distances(num_nodes_, std::numeric_limits<Cost>::infinity());
  PredecessorVector predecessors(num_nodes_, NO_NODE);

  if (settle_order)
  {
    settle_order->clear();
    settle_order->reserve(num_nodes_);
  }

  using QueueNode = std::pair<Cost, NodeId>;
  std::priority_queue<QueueNode, std::vector<QueueNode>, std::greater<>> open_set;

  for (const Source& s : sources)
  {
    distances[s.node] = s.cost;
    open_set.push({s.cost, s.node});
  }

  Counter counter("Dijkstra", num_nodes_);
  std::size_t settled = 0;

  while (!open_set.empty())
  {
    // No structured bindings: capturing them in a lambda is C++20 but unsupported by older clang (macOS oldrel)
    const Cost current_cost = open_set.top().first;
    const NodeId current = open_set.top().second;
    open_set.pop();

    if (current_cost > distances[current]) continue;

    counter.tick_at(settled++);
    if (counter.cancelled()) break;

    // A node is pushed only on a strict improvement, so this point is reached once per node
    if (settle_order) settle_order->push_back(current);

    for_each_edge(current, [&](const Node& neighbor)
    {
      Cost new_cost = current_cost + neighbor.cost;
      if (new_cost < distances[neighbor.destination])
      {
        distances[neighbor.destination] = new_cost;
        predecessors[neighbor.destination] = current;
        open_set.push({new_cost, neighbor.destination});
      }
    });
  }

  counter.throw_if_cancelled();

  return {std::move(distances), std::move(predecessors)};
}

void Graph::shortest_paths_from_sources(std::vector<double>& distances, NodeIDs& closest_seed) const
{
  const double INF = std::numeric_limits<double>::infinity();

  distances.assign(num_nodes_, INF);
  closest_seed.assign(num_nodes_, NO_NODE);

  using PQElement = std::pair<double, NodeId>;
  std::priority_queue<PQElement, std::vector<PQElement>, std::greater<>> pq;

  for (const Source& s : sources)
  {
    distances[s.node] = s.cost;
    closest_seed[s.node] = s.seed;
    pq.push({distances[s.node], s.node});
  }

  Counter counter("Dijkstra", num_nodes_);
  std::size_t settled = 0;

  while (!pq.empty())
  {
    const double dist_u = pq.top().first;
    const NodeId u = pq.top().second;
    pq.pop();

    if (dist_u > distances[u]) continue;

    counter.tick_at(settled++);
    if (counter.cancelled()) break;

    for_each_edge(u, [&](const Node& e)
    {
      NodeId v = e.destination;
      double new_dist = dist_u + e.cost;
      if (new_dist < distances[v])
      {
        distances[v] = new_dist;
        closest_seed[v] = closest_seed[u];
        pq.push({new_dist, v});
      }
    });
  }

  counter.throw_if_cancelled();
}

size_t Graph::mem() const
{
  size_t total = sizeof(*this);
  total += num_nodes_ * stride_ * sizeof(Node);
  total += extra_offset_.capacity() * sizeof(uint32_t);
  total += extra_edges_.capacity() * sizeof(Node);
  total += sources.capacity() * sizeof(Source);
  total += target_anchors.capacity() * sizeof(NodeId);
  return total;
}

size_t Graph::cache_mem(const GraphCache& cache)
{
  const auto& [distances, predecessors] = cache;
  size_t total = sizeof(distances) + distances.capacity() * sizeof(Cost);
  total += predecessors.size() * sizeof(NodeId);
  return total;
}

std::string Graph::format_bytes(size_t bytes)
{
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2);

  if (bytes < 1024)
    oss << bytes << " B";
  else if (bytes < 1024 * 1024)
    oss << (bytes / 1024.0) << " KB";
  else if (bytes < 1024 * 1024 * 1024)
    oss << (bytes / (1024.0 * 1024.0)) << " MB";
  else
    oss << (bytes / (1024.0 * 1024.0 * 1024.0)) << " GB";

  return oss.str();
}

}

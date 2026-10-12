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

#include "Graph.h"

namespace arbor::segment {

Graph::Graph(size_t num_nodes, size_t stride) :
  num_nodes_(num_nodes),
  stride_(stride),
  slots_(num_nodes * stride, Node{NO_NODE, 0.0f}),
  has_extra_(num_nodes, false),
  extra_offset_(1, 0)
{
}

void Graph::clear_row(NodeId u)
{
  Node* r = row(u);
  std::fill(r, r + stride_, Node{NO_NODE, 0.0f});
}

void Graph::add_extra_edges(const std::vector<Edge>& edges)
{
  if (edges.empty()) return;

  // Rebuild the side table from the existing and the new edges. This is called at most a
  // couple of times on a small number of edges.
  std::vector<Edge> all;
  all.reserve(extra_edges_.size() + edges.size());
  for (size_t i = 0; i < extra_src_.size(); ++i)
  {
    for (size_t e = extra_offset_[i]; e < extra_offset_[i + 1]; ++e)
      all.push_back({extra_src_[i], extra_edges_[e]});
  }
  all.insert(all.end(), edges.begin(), edges.end());

  // Stable: keeps the insertion order of the edges of a given node
  std::stable_sort(all.begin(), all.end(), [](const Edge& a, const Edge& b) { return a.source < b.source; });

  extra_src_.clear();
  extra_offset_.assign(1, 0);
  extra_edges_.clear();
  extra_edges_.reserve(all.size());

  for (const Edge& e : all)
  {
    if (extra_src_.empty() || extra_src_.back() != e.source)
    {
      if (!extra_src_.empty()) extra_offset_.push_back(extra_edges_.size());
      extra_src_.push_back(e.source);
      has_extra_[e.source] = true;
    }
    extra_edges_.push_back(e.node);
  }
  extra_offset_.push_back(extra_edges_.size());
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

  while (!open_set.empty())
  {
    auto [current_cost, current] = open_set.top();
    open_set.pop();

    if (current_cost > distances[current]) continue;

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

  while (!pq.empty())
  {
    auto [dist_u, u] = pq.top();
    pq.pop();

    if (dist_u > distances[u]) continue;

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
}

size_t Graph::mem() const
{
  size_t total = sizeof(*this);
  total += slots_.capacity() * sizeof(Node);
  total += has_extra_.capacity() / 8;
  total += extra_src_.capacity() * sizeof(NodeId);
  total += extra_offset_.capacity() * sizeof(size_t);
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

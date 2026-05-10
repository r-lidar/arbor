/**
 * @file Graph.h
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

#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <cstddef>
#include <sstream>
#include <iomanip>

namespace arbor::segment {

class Graph
{
public:
  using NodeId = int;
  using Cost = float;
  using DistanceVector = std::vector<Cost>;
  using PredecessorVector = std::vector<NodeId>;
  using Path = std::vector<NodeId>;
  using NodeIDs = std::vector<NodeId>;
  using GraphCache = std::pair<DistanceVector, PredecessorVector>;

  struct Node
  {
    NodeId destination;
    Cost cost;
  };

  using AdjacencyList = std::vector<std::vector<Node>>;

public:
  AdjacencyList adjacency_list;

  Graph() = default;
  void ensure_size(size_t n);
  void reserve_edges(size_t capacity);
  void add_edge(NodeId source, NodeId destination, Cost cost);
  void clear_adjacency_list();
  std::pair<DistanceVector, PredecessorVector> compute_distances(NodeId start) const;
  std::pair<Path, Cost> findPath(NodeId start, NodeId goal, const std::pair<DistanceVector, PredecessorVector>& precomputed_data) const;
  void shortest_paths_from_node(const NodeIDs& ground_nodes, std::vector<double>& distances, NodeIDs& closest_ground) const;

  static size_t cache_mem(const Graph::GraphCache& cache)
  {
    const auto& [distances, predecessors] = cache;
    size_t total = sizeof(distances) + distances.capacity() * sizeof(Cost);
    total += predecessors.size() * sizeof(NodeId);
    return total;
  }

  // Calculate total memory size of a Graph in bytes
  size_t mem()
  {
    size_t total = sizeof(this);
    for (const auto& neighbors : adjacency_list)
    {
      total += sizeof(neighbors); // vector overhead per node
      total += neighbors.capacity() * sizeof(Node);
    }

    return total;
  }

  // Format bytes into human-readable string (KB, MB, GB)
  static std::string format_bytes(size_t bytes)
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
};

}

#endif // GRAPH_H

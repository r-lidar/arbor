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
#include <unordered_map>
#include <cstddef>

namespace arbor::segment {

class Graph
{
public:
  using NodeId = int;
  using Cost = float;
  using DistanceVector = std::vector<Cost>;
  using PredecessorMap = std::unordered_map<NodeId, NodeId>;
  using Path = std::vector<NodeId>;
  using NodeIDs = std::vector<NodeId>;
  using GraphCache = std::pair<DistanceVector, PredecessorMap>;

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
  void add_edge(NodeId source, NodeId destination, Cost cost);
  std::pair<DistanceVector, PredecessorMap> compute_distances(NodeId start) const;
  std::pair<Path, Cost> findPath(NodeId start, NodeId goal, const std::pair<DistanceVector, PredecessorMap>& precomputed_data) const;
  void shortest_paths_from_node(const NodeIDs& ground_nodes, std::vector<double>& distances, NodeIDs& closest_ground) const;
};

}

#endif // GRAPH_H

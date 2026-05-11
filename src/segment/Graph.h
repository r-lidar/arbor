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

// The graph uses a two-phase Compressed Sparse Row (CSR) representation to
// minimise memory use and improve cache locality for large graphs (25M+ nodes).
//
// Build phase  : edges are accumulated in a flat COO buffer (coo_src_ / coo_edges_).
// Query phase  : on the first call to compute_distances() or
//                shortest_paths_from_node(), finalize() converts COO → CSR and
//                frees the COO buffers.
// Overflow     : edges added *after* CSR finalisation (rare fallback in dist2root)
//                are stored in a small hash-map and merged into the traversal
//                with O(1) per-node lookup.

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

public:
  Graph() = default;
  void ensure_size(size_t n);
  void add_edge(NodeId source, NodeId destination, Cost cost);
  std::pair<DistanceVector, PredecessorMap> compute_distances(NodeId start) const;
  std::pair<Path, Cost> findPath(NodeId start, NodeId goal, const std::pair<DistanceVector, PredecessorMap>& precomputed_data) const;
  void shortest_paths_from_node(const NodeIDs& ground_nodes, std::vector<double>& distances, NodeIDs& closest_ground) const;

private:
  size_t num_nodes_ = 0;

  // COO build buffers — populated by add_edge(), freed by finalize().
  // Marked mutable to support the lazy-initialisation pattern: compute_distances()
  // and shortest_paths_from_node() are logically const (they do not change the
  // observable graph) but must trigger finalize() on the first call to convert
  // the COO buffers to CSR format.  This is the standard "mutable cache" idiom.
  mutable std::vector<NodeId> coo_src_;
  mutable std::vector<Node>   coo_edges_;

  // CSR storage — built lazily by finalize().
  // row_ptr_[i] .. row_ptr_[i+1] is the half-open range of adj_data_ for node i.
  mutable std::vector<size_t> row_ptr_;
  mutable std::vector<Node>   adj_data_;
  mutable bool csr_valid_ = false;

  // Overflow: edges added after CSR has been finalised (rare).
  // Kept as a small adjacency map so Dijkstra can look them up in O(1).
  std::unordered_map<NodeId, std::vector<Node>> overflow_;

  void finalize() const; // converts COO → CSR; no-op if already done
};

}

#endif // GRAPH_H

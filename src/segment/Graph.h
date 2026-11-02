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

#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <cstddef>

namespace arbor::segment {

// Directed graph over the core points only.
//
// The core layer is a kNN graph: each node has at most `stride` out-edges. They are stored in
// one flat array of `stride` slots per node (unused slots hold NO_NODE). The extra edges added
// by the reachability repair are kept in a side table in CSR layout, indexed by source node.
//
// The seed and master layers of the original layered graph are collapsed into `sources`: the
// core nodes reached by a seed, with their initial cost and seed index (multi-source Dijkstra).
// The target layer is collapsed into `target_anchors`: each target is a leaf hanging off its
// nearest core node with a zero cost, so its distance and path are those of its anchor.
class Graph
{
public:
  using NodeId = unsigned int;
  using Cost = float;
  using DistanceVector = std::vector<Cost>;
  using PredecessorVector = std::vector<NodeId>;
  using NodeIDs = std::vector<NodeId>;
  using GraphCache = std::pair<DistanceVector, PredecessorVector>;

  static constexpr NodeId NO_NODE = std::numeric_limits<NodeId>::max();

  struct Node
  {
    NodeId destination;
    Cost cost;
  };

  struct Edge
  {
    NodeId source;
    Node node;
  };

  struct Source
  {
    NodeId node;  // core node reached by a seed
    Cost cost;    // cost from the seed to this node
    NodeId seed;  // index of the seed
  };

public:
  // Collapsed seed and master layers
  std::vector<Source> sources;

  // Collapsed target layer: nearest core node of each target, NO_NODE if the target was dropped
  NodeIDs target_anchors;

  // The slots are allocated but not initialized: every slot of every row must be written by the
  // builder (edges first, then NO_NODE in the unused slots).
  Graph(size_t num_nodes, size_t stride);

  size_t size() const { return num_nodes_; }
  size_t stride() const { return stride_; }

  // Fixed-size row of out-edges of node u, filled from the start, unused slots hold NO_NODE
  Node* row(NodeId u) { return &slots_[static_cast<size_t>(u) * stride_]; }
  void clear_row(NodeId u);

  // Adds edges to the side table and releases the input. Edges of a node are visited after its
  // row, in insertion order.
  void add_extra_edges(std::vector<Edge> edges);

  template <typename F>
  void for_each_edge(NodeId u, F&& f) const
  {
    const Node* r = &slots_[static_cast<size_t>(u) * stride_];
    for (size_t j = 0; j < stride_ && r[j].destination != NO_NODE; ++j)
      f(r[j]);

    if (extra_offset_.empty()) return;

    for (uint32_t e = extra_offset_[u]; e < extra_offset_[u + 1]; ++e)
      f(extra_edges_[e]);
  }

  // Dijkstra from all the sources. Roots (source nodes not improved by another path) have no
  // predecessor. If settle_order is provided, it receives the reached nodes in the order Dijkstra
  // settled them. This is a topological order of the shortest-path tree: every node comes after
  // its predecessor.
  GraphCache compute_distances(NodeIDs* settle_order = nullptr) const;

  // Dijkstra from all the sources. closest_seed receives for each node the index of the seed it
  // is the closest to, or NO_NODE if unreachable.
  void shortest_paths_from_sources(std::vector<double>& distances, NodeIDs& closest_seed) const;

  size_t mem() const;
  static size_t cache_mem(const GraphCache& cache);
  static std::string format_bytes(size_t bytes);

private:
  size_t num_nodes_;
  size_t stride_;
  std::unique_ptr<Node[]> slots_;

  // Side table of extra edges in CSR layout: edges of node u are
  // extra_edges_[extra_offset_[u] .. extra_offset_[u+1]). Both are empty until edges are added.
  std::vector<uint32_t> extra_offset_;
  std::vector<Node> extra_edges_;
};

}

#endif // GRAPH_H

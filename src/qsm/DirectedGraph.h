/**
 * @file DirectedGraph.h
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

#ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <stdexcept>

using NodeID = int;
using EdgeID = int;

// Generic directed graph with typed node and edge data.
// NodeID and EdgeID are integers, auto-incremented from 1.
template<typename NodeData, typename EdgeData>
class DirectedGraph
{
public:
  struct EdgeInfo
  {
    NodeID  source;
    NodeID  target;
    EdgeData data;
  };

  // ---- Node operations ----

  NodeID add_node(const NodeData& data)
  {
    NodeID id = next_node_id_++;
    nodes_[id] = data;
    outgoing_[id] = {};
    incoming_[id] = {};
    return id;
  }

  bool has_node(NodeID id) const { return nodes_.count(id) > 0; }
  size_t node_count() const { return nodes_.size(); }

  NodeData& node(NodeID id) { return nodes_.at(id); }
  const NodeData& node(NodeID id) const { return nodes_.at(id); }

  void remove_node(NodeID id)
  {
    nodes_.erase(id);
    outgoing_.erase(id);
    incoming_.erase(id);
  }

  // ---- Edge operations ----

  EdgeID add_edge(NodeID source, NodeID target, const EdgeData& data)
  {
    EdgeID id = next_edge_id_++;
    edges_[id] = {source, target, data};
    outgoing_[source].push_back(id);
    incoming_[target].push_back(id);
    return id;
  }

  bool has_edge(EdgeID id) const { return edges_.count(id) > 0; }
  size_t edge_count() const { return edges_.size(); }

  EdgeInfo& edge(EdgeID id) { return edges_.at(id); }
  const EdgeInfo& edge(EdgeID id) const { return edges_.at(id); }

  EdgeData& edge_data(EdgeID id) { return edges_.at(id).data; }
  const EdgeData& edge_data(EdgeID id) const { return edges_.at(id).data; }

  void remove_edge(EdgeID id)
  {
    auto it = edges_.find(id);
    if (it == edges_.end()) return;

    NodeID src = it->second.source;
    NodeID tgt = it->second.target;

    auto& out = outgoing_[src];
    out.erase(std::remove(out.begin(), out.end(), id), out.end());

    auto& inc = incoming_[tgt];
    inc.erase(std::remove(inc.begin(), inc.end(), id), inc.end());

    edges_.erase(it);
  }

  // ---- Graph traversal ----

  const std::vector<EdgeID>& outgoing_edges(NodeID id) const
  {
    static const std::vector<EdgeID> empty{};
    auto it = outgoing_.find(id);
    return it != outgoing_.end() ? it->second : empty;
  }

  const std::vector<EdgeID>& incoming_edges(NodeID id) const
  {
    static const std::vector<EdgeID> empty{};
    auto it = incoming_.find(id);
    return it != incoming_.end() ? it->second : empty;
  }

  std::vector<NodeID> children(NodeID id) const
  {
    std::vector<NodeID> result;
    for (EdgeID eid : outgoing_edges(id))
      result.push_back(edges_.at(eid).target);
    return result;
  }

  // Returns the source node of the first incoming edge (-1 if none).
  NodeID parent_node(NodeID id) const
  {
    const auto& inc = incoming_edges(id);
    if (inc.empty()) return -1;
    return edges_.at(inc[0]).source;
  }

  // ---- Storage access ----

  const std::unordered_map<NodeID, NodeData>& nodes() const { return nodes_; }
  std::unordered_map<NodeID, NodeData>&       nodes()       { return nodes_; }

  const std::unordered_map<EdgeID, EdgeInfo>& edges() const { return edges_; }
  std::unordered_map<EdgeID, EdgeInfo>&       edges()       { return edges_; }

  void clear()
  {
    nodes_.clear();
    edges_.clear();
    outgoing_.clear();
    incoming_.clear();
    next_node_id_ = 1;
    next_edge_id_ = 1;
  }

private:
  std::unordered_map<NodeID, NodeData>              nodes_;
  std::unordered_map<EdgeID, EdgeInfo>              edges_;
  std::unordered_map<NodeID, std::vector<EdgeID>>   outgoing_;
  std::unordered_map<NodeID, std::vector<EdgeID>>   incoming_;

  NodeID next_node_id_ = 1;
  EdgeID next_edge_id_ = 1;
};

#endif // DIRECTED_GRAPH_H

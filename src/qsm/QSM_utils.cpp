/**
 * @file QSM_utils.cpp
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

#include "QSM.h"

#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <string>

namespace arbor::qsm {

void QSM::validate() const
{
  // ------------------------------------------------------------------ //
  // 1. Single-root check
  // ------------------------------------------------------------------ //
  NodeID root      = -1;
  int    root_count = 0;

  for (const auto& [id, _] : nodes())
  {
    if (incoming_edges(id).empty())
    {
      ++root_count;
      root = id;
    }
  }

  if (root_count == 0)
    throw std::runtime_error("QSM::validate: graph has no root node.");

  if (root_count > 1)
    throw std::runtime_error("QSM::validate: graph has " + std::to_string(root_count) + " root nodes.");

  // ------------------------------------------------------------------ //
  // 2. Unset-radius check
  // ------------------------------------------------------------------ //
  for (const auto& [eid, info] : edges())
  {
    if (info.data.radius == RADIUS_UNSET)
      throw std::runtime_error("QSM::validate: edge " + std::to_string(eid) + " has NA radius.");

    if (std::isnan(info.data.radius))
      throw std::runtime_error("QSM::validate: edge " + std::to_string(eid) + " has NaN radius");
  }

  // ------------------------------------------------------------------ //
  // 3. Connectivity check — BFS from root
  // ------------------------------------------------------------------ //
  std::unordered_set<NodeID> visited;
  std::vector<NodeID>        stack = { root };

  while (!stack.empty())
  {
    NodeID cur = stack.back();
    stack.pop_back();
    if (!visited.insert(cur).second) continue;

    for (EdgeID eid : outgoing_edges(cur))
      stack.push_back(edge(eid).target);
  }

  for (const auto& [id, _] : nodes())
  {
    if (visited.find(id) == visited.end())
      throw std::runtime_error("QSM::validate: node " + std::to_string(id) + " is not reachable from the root.");
  }
}

NodeID QSM::find_root_node() const
{
  for (const auto& [nid, _] : nodes())
  {
    if (incoming_edges(nid).empty())
      return nid;
  }
  return static_cast<NodeID>(-1);
}

// Returns a new QSM containing only the edges (and their endpoint nodes)
// whose axis_id == 1 (i.e. the main trunk axis).
// Node IDs are remapped; distance_to_root values are preserved as-is.
QSM QSM::stem() const
{
  QSM result;

  // Map from original NodeID -> new NodeID in the result graph.
  std::unordered_map<NodeID, NodeID> node_map;

  auto get_or_add_node = [&](NodeID orig) -> NodeID
  {
    auto it = node_map.find(orig);
    if (it != node_map.end()) return it->second;
    NodeID new_id = result.add_node(node(orig));
    node_map[orig] = new_id;
    return new_id;
  };

  for (const auto& [eid, einfo] : edges())
  {
    if (einfo.data.branch_order == 1)
    {
      NodeID new_src = get_or_add_node(einfo.source);
      NodeID new_tgt = get_or_add_node(einfo.target);
      result.add_edge(new_src, new_tgt, einfo.data);
    }
  }

  return result;
}

QSM QSM::merchantable(double min_radius, double min_axis_length) const
{
  QSM result;
  std::unordered_set<EdgeID> edges_to_keep;

  // STEP 1: Radius-based Pruning from leaves to root. Stop at first edge with big enought radius
  std::vector<NodeID> leaves;
  for (const auto& [node_id, _] : nodes())
  {
    if (outgoing_edges(node_id).empty())
    {
      leaves.push_back(node_id);
    }
  }

  for (NodeID leaf_id : leaves)
  {
    NodeID current = leaf_id;
    while (current != -1)
    {
      const auto& inc_edges = incoming_edges(current);
      if (inc_edges.empty()) break;

      EdgeID parent_edge_id = inc_edges[0];
      const auto& edge_info = edge(parent_edge_id);
      if (edge_info.data.radius > min_radius)
      {
        // Keep this edge and traverse back to root
        NodeID node_on_path = current;
        while (node_on_path != -1)
        {
          const auto& path_inc_edges = incoming_edges(node_on_path);
          if (path_inc_edges.empty()) break;

          EdgeID path_edge_id = path_inc_edges[0];
          edges_to_keep.insert(path_edge_id);
          node_on_path = edge(path_edge_id).source;
        }
        break;
      }
      current = edge_info.source;
    }
  }

  // STEP 2: Axis Length Filtering
  // Calculate total length for each remaining axis
  std::unordered_map<int, double> axis_total_lengths;
  for (EdgeID eid : edges_to_keep)
  {
    const auto& e_info = edge(eid);
    double len = e_info.data.length(node(e_info.source), node(e_info.target));
    axis_total_lengths[e_info.data.axis_id] += len;
  }

  // Filter out edges belonging to axes that are too short
  for (auto it = edges_to_keep.begin(); it != edges_to_keep.end(); )
  {
    const auto& e_info = edge(*it);
    if (axis_total_lengths[e_info.data.axis_id] < min_axis_length)
    {
      it = edges_to_keep.erase(it);
    }
    else
    {
      ++it;
    }
  }

  // STEP 3: Graph Reconstruction
  std::unordered_set<NodeID> nodes_to_keep;
  for (EdgeID edge_id : edges_to_keep)
  {
    const auto& edge_info = edge(edge_id);
    nodes_to_keep.insert(edge_info.source);
    nodes_to_keep.insert(edge_info.target);
  }

  std::unordered_map<NodeID, NodeID> old_to_new;
  for (NodeID old_id : nodes_to_keep)
  {
    old_to_new[old_id] = result.add_node(node(old_id));
  }

  for (EdgeID edge_id : edges_to_keep)
  {
    const auto& edge_info = edge(edge_id);
    result.add_edge(old_to_new[edge_info.source], old_to_new[edge_info.target], edge_info.data);
  }

  return result;
}

double QSM::dbh(double d, double* xyz, double* n) const
{
  // Collect trunk edges sorted by distance_to_root.
  struct TrunkEdge
  {
    double distance_to_root;
    double radius;
    NodeID source;
    NodeID target;
  };

  std::vector<TrunkEdge> trunk_edges;
  trunk_edges.reserve(64);

  // Collect trunk edges
  for (const auto& [eid, einfo] : edges())
  {
    if (einfo.data.axis_id != 1) continue;
    if (einfo.data.distance_to_root == DISTANCE_TO_ROOT_UNSET) throw std::runtime_error("'Distance to root' attribute not populated");
    if (einfo.data.distance_to_root == RADIUS_UNSET) throw std::runtime_error("'Radius' attribute not populated");
    trunk_edges.push_back({einfo.data.distance_to_root, einfo.data.radius, einfo.source, einfo.target});
  }

  if (trunk_edges.empty()) throw std::runtime_error("Internal error: no axis_id = 1 in this QSM");

  // Sort by distance to root
  std::sort(trunk_edges.begin(), trunk_edges.end(), [](const TrunkEdge& a, const TrunkEdge& b) { return a.distance_to_root < b.distance_to_root; });

  // Helper: fill xyz and n from a single edge (clamp case).
  auto fill_from_single = [&](const TrunkEdge& e)
  {
    const QSMNode& src = node(e.source);
    const QSMNode& tgt = node(e.target);

    // Position: midpoint of the edge.
    if (xyz)
    {
      xyz[0] = 0.5 * (src.x + tgt.x);
      xyz[1] = 0.5 * (src.y + tgt.y);
      xyz[2] = 0.5 * (src.z + tgt.z);
    }

    // Normal: direction of the edge, normalised.
    if (n)
    {
      double dx = tgt.x - src.x;
      double dy = tgt.y - src.y;
      double dz = tgt.z - src.z;
      double len = std::sqrt(dx*dx + dy*dy + dz*dz);
      if (len < 1e-12) len = 1e-12;
      n[0] = dx / len;
      n[1] = dy / len;
      n[2] = dz / len;
    }
  };

  // Clamp to the range covered by trunk edges.
  if (d <= trunk_edges.front().distance_to_root)
  {
    fill_from_single(trunk_edges.front());
    return trunk_edges.front().radius;
  }
  if (d >= trunk_edges.back().distance_to_root)
  {
    fill_from_single(trunk_edges.back());
    return trunk_edges.back().radius;
  }

  // Binary search for the first edge whose distance_to_root > d.
  auto it = std::upper_bound(trunk_edges.begin(), trunk_edges.end(), d, [](double val, const TrunkEdge& e)
  {
    return val < e.distance_to_root;
  });

  const TrunkEdge& after  = *it;
  const TrunkEdge& before = *std::prev(it);

  double span = after.distance_to_root - before.distance_to_root;
  double t    = (span < 1e-12) ? 0.0 : (d - before.distance_to_root) / span;

  // Interpolated radius.
  double radius = before.radius + t * (after.radius - before.radius);

  // Interpolated position:
  if (xyz)
  {
    const QSMNode& p0 = node(before.target);
    const QSMNode& p1 = node(after.source);
    xyz[0] = p0.x + t * (p1.x - p0.x);
    xyz[1] = p0.y + t * (p1.y - p0.y);
    xyz[2] = p0.z + t * (p1.z - p0.z);
  }

  // Normal: direction vector from before.source -> after.target, normalised.
  if (n)
  {
    const QSMNode& p0 = node(before.source);
    const QSMNode& p1 = node(after.target);
    double dx = p1.x - p0.x;
    double dy = p1.y - p0.y;
    double dz = p1.z - p0.z;
    double len = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (len < 1e-12) len = 1e-12;
    n[0] = dx / len;
    n[1] = dy / len;
    n[2] = dz / len;
  }

  return radius*2;
}



}

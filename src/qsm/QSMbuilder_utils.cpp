/**
 * @file QSMbuilder_utils.cpp
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

#include "QSMbuilder.h"

namespace arbor::qsm {

void QSMbuilder::shift(double tx, double ty, double tz)
{
  ServiceLocator::logger()("Shift back to geographic coordinates");

  // Shift all node positions (a single update per node covers all incident edges)
  for (auto& [nid, ndata] : graph.nodes())
  {
    ndata.x += tx;
    ndata.y += ty;
    ndata.z += tz;
  }
}

int QSMbuilder::count_nodes_connected_to_root() const
{
  // Iterate through the public nodes map
  for (const auto& [id, data] : graph.nodes())
  {
    // The root is defined as a node with zero incoming edges
    if (graph.incoming_edges(id).empty())
    {
      // Return the count of its outgoing edges
      return graph.outgoing_edges(id).size();
    }
  }

  return 0;
}

int QSMbuilder::find_root_edge() const
{
  for (const auto& [eid, einfo] : graph.edges())
  {
    if (graph.incoming_edges(einfo.source).empty())
      return eid;
  }
  return -1;
}

void QSM::repair_multiroot()
{
  const std::vector<NodeID> roots = this->roots();
  if (roots.size() <= 1) return;

  struct Score
  {
    double volume = 0.0;
    double length = 0.0;
    size_t nodes  = 0;
  };

  auto better = [](const Score& a, const Score& b)
  {
    if (a.volume != b.volume) return a.volume > b.volume;
    if (a.length != b.length) return a.length > b.length;
    return a.nodes > b.nodes;
  };

  // Iterative traversal of everything reachable from `root`.
  auto explore = [this](NodeID root, std::unordered_set<NodeID>& reached) -> Score
  {
    Score s;
    std::vector<NodeID> stack{root};
    reached.insert(root);

    while (!stack.empty())
    {
      NodeID n = stack.back();
      stack.pop_back();
      s.nodes++;

      for (EdgeID eid : outgoing_edges(n))
      {
        const auto& e = edge(eid);
        const QSMNode& src = node(e.source);
        const QSMNode& tgt = node(e.target);
        s.volume += e.data.volume(src, tgt);
        s.length += e.data.length(src, tgt);

        if (reached.insert(e.target).second)
          stack.push_back(e.target);
      }
    }
    return s;
  };

  // Find the best component
  Score best;
  NodeID best_root = -1;
  std::unordered_set<NodeID> kept;

  for (NodeID r : roots)
  {
    std::unordered_set<NodeID> reached;
    Score s = explore(r, reached);

    if (best_root == -1 || better(s, best))
    {
      best = s;
      best_root = r;
      kept = std::move(reached);
    }
  }

  const size_t edges_before = edge_count();

  std::vector<NodeID> nodes_to_remove;
  for (const auto& kv : nodes())
  {
    if (!kept.count(kv.first))
      nodes_to_remove.push_back(kv.first);
  }

  for (NodeID nid : nodes_to_remove)
    remove_node(nid);

  const size_t edges_removed = edges_before - edge_count();

  messages.push_back("[E2] [Topology issue] The tree model had disconnected parts. Auto-repaired by keeping only the main part");

  return;
}

} // namespace arbor::qsm

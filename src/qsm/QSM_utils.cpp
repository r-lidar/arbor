#include "QSM.h"

#include <unordered_set>

namespace arbor::qsm {

NodeID QSM::find_root_node() const
{
  for (const auto& [nid, _] : nodes())
    if (incoming_edges(nid).empty())
      return nid;
    return -1;
}

// Returns a new QSM containing only the edges (and their endpoint nodes)
// whose axis_ID == 1 (i.e. the main trunk axis).
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

QSM QSM::merchantable(double min_radius) const
{
  QSM result;

  // Track which edges should be kept
  std::unordered_set<EdgeID> edges_to_keep;

  // Find all leaf nodes (nodes with no outgoing edges)
  std::vector<NodeID> leaves;
  for (const auto& [node_id, _] : nodes())
  {
    if (outgoing_edges(node_id).empty())
    {
      leaves.push_back(node_id);
    }
  }

  // For each leaf, traverse back to root
  for (NodeID leaf_id : leaves)
  {
    NodeID current = leaf_id;

    // Traverse from leaf to root
    while (current != -1)
    {
      const auto& inc_edges = incoming_edges(current);

      // If no incoming edges, we've reached the root
      if (inc_edges.empty())
      {
        break;
      }

      // Get the parent edge (first incoming edge)
      EdgeID parent_edge_id = inc_edges[0];
      const auto& edge_info = edge(parent_edge_id);

      // Check if this edge meets the radius threshold
      if (edge_info.data.radius >= min_radius)
      {
        // This edge is valid - keep it and all edges back to root
        NodeID node_on_path = current;

        while (node_on_path != -1)
        {
          const auto& path_inc_edges = incoming_edges(node_on_path);

          if (path_inc_edges.empty())
          {
            break;
          }

          EdgeID path_edge_id = path_inc_edges[0];
          edges_to_keep.insert(path_edge_id);

          const auto& path_edge_info = edge(path_edge_id);
          node_on_path = path_edge_info.source;
        }

        break; // Done with this leaf's path
      }

      // Edge is too thin - continue moving toward root
      current = edge_info.source;
    }
  }

  // Collect all nodes that are endpoints of kept edges
  std::unordered_set<NodeID> nodes_to_keep;
  for (EdgeID edge_id : edges_to_keep)
  {
    const auto& edge_info = edge(edge_id);
    nodes_to_keep.insert(edge_info.source);
    nodes_to_keep.insert(edge_info.target);
  }

  // Build the new graph
  std::unordered_map<NodeID, NodeID> old_to_new;

  // Add nodes
  for (NodeID old_id : nodes_to_keep)
  {
    NodeID new_id = result.add_node(node(old_id));
    old_to_new[old_id] = new_id;
  }

  // Add edges
  for (EdgeID edge_id : edges_to_keep)
  {
    const auto& edge_info = edge(edge_id);
    NodeID new_source = old_to_new[edge_info.source];
    NodeID new_target = old_to_new[edge_info.target];
    result.add_edge(new_source, new_target, edge_info.data);
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
    if (einfo.data.axis_ID != 1) continue;
    if (einfo.data.distance_to_root == DISTANCE_TO_ROOT_UNSET) throw std::runtime_error("'Distance to root' attribute not populated");
    if (einfo.data.distance_to_root == RADIUS_UNSET) throw std::runtime_error("'Radius' attribute not populated");
    trunk_edges.push_back({einfo.data.distance_to_root, einfo.data.radius, einfo.source, einfo.target});
  }

  if (trunk_edges.empty()) throw std::runtime_error("Internal error: no axis_ID = 1 in this QSM");

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

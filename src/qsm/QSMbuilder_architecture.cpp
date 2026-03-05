#include "QSMbuilder.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace arbor::qsm {

void QSMbuilder::compute_topology()
{
  logger("Connecting nodes");

  // For each edge, set parent_ID based on graph structure:
  //   - If the edge's source node has no incoming edges -> parent_ID = 0 (root)
  //   - Otherwise -> parent_ID = cyl_ID of the (first) incoming edge to source
  for (auto& [eid, einfo] : graph.edges())
  {
    const auto& inc = graph.incoming_edges(einfo.source);
    if (inc.empty())
    {
      einfo.data.parent_ID = 0;
    }
    else
    {
      einfo.data.parent_ID = graph.edge_data(inc[0]).cyl_ID;
    }
  }
}

void QSMbuilder::compute_architecture(int /*root_id*/, bool use_volume)
{
  logger("Computing architecture");

  // Reset cached values on every edge
  for (auto& [eid, einfo] : graph.edges())
  {
    einfo.data.subtree_length   = SUBTREE_LENGTH_UNSET;
    einfo.data.subtree_max_endZ = SUBTREE_MAXZ_UNSET;
    einfo.data.subtree_volume   = SUBTREE_VOLUME_UNSET;
    einfo.data.axis_ID          = 0;
    einfo.data.branch_order     = 0;
  }

  // Find the root edge (source has no incoming edges)
  int root_eid = find_root_edge();
  if (root_eid < 0) return;

  compute_subtree_length(root_eid);

  if (use_volume)
    compute_subtree_volume(root_eid);
  else
    compute_subtree_max_z(root_eid);

  int next_axis_id = 2;
  assign_subtree_ids(root_eid, 1, 1, next_axis_id, use_volume);
}

double QSMbuilder::compute_subtree_length(int edge_id)
{
  auto& ed = graph.edge_data(edge_id);

  if (ed.subtree_length >= 0)
    return ed.subtree_length;

  const QSMNode& src = graph.node(graph.edge(edge_id).source);
  const QSMNode& tgt = graph.node(graph.edge(edge_id).target);

  const auto& child_eids = graph.outgoing_edges(graph.edge(edge_id).target);
  if (child_eids.empty())
  {
    ed.subtree_length = 0.0;
    return 0.0;
  }

  double max_len = 0.0;
  for (int child_eid : child_eids)
  {
    const auto& child_einfo = graph.edge(child_eid);
    double child_len = child_einfo.data.length(
      graph.node(child_einfo.source),
      graph.node(child_einfo.target));
    double candidate = compute_subtree_length(child_eid) + child_len;
    max_len = std::max(max_len, candidate);
  }

  ed.subtree_length = max_len;
  return max_len;
}

double QSMbuilder::compute_subtree_max_z(int edge_id)
{
  auto& ed = graph.edge_data(edge_id);

  if (ed.subtree_max_endZ > SUBTREE_MAXZ_UNSET)
    return ed.subtree_max_endZ;

  double maxz = graph.node(graph.edge(edge_id).target).z;

  for (int child_eid : graph.outgoing_edges(graph.edge(edge_id).target))
  {
    double child_maxz = compute_subtree_max_z(child_eid);
    maxz = std::max(maxz, child_maxz);
  }

  ed.subtree_max_endZ = maxz;
  return maxz;
}

double QSMbuilder::compute_subtree_volume(int edge_id)
{
  auto& ed = graph.edge_data(edge_id);

  if (ed.subtree_volume >= 0)
    return ed.subtree_volume;

  const auto& child_eids = graph.outgoing_edges(graph.edge(edge_id).target);
  if (child_eids.empty())
  {
    ed.subtree_volume = 0.0;
    return 0.0;
  }

  double max_v = 0.0;
  for (int child_eid : child_eids)
  {
    const auto& child_einfo = graph.edge(child_eid);
    double child_vol = child_einfo.data.volume(
      graph.node(child_einfo.source),
      graph.node(child_einfo.target));
    double candidate = compute_subtree_volume(child_eid) + child_vol;
    max_v = std::max(max_v, candidate);
  }

  ed.subtree_volume = max_v;
  return max_v;
}

void QSMbuilder::assign_subtree_ids(int edge_id, int current_axis_id, int current_branch_order, int& next_axis_id, bool use_volume)
{
  auto& ed = graph.edge_data(edge_id);
  ed.axis_ID      = current_axis_id;
  ed.branch_order = current_branch_order;

  const auto& child_eids = graph.outgoing_edges(graph.edge(edge_id).target);
  if (child_eids.empty()) return;

  int main_child = -1;
  double best_score     = -1e300;
  double best_secondary = -1e300;

  for (int child_eid : child_eids)
  {
    const auto& child_einfo = graph.edge(child_eid);
    const QSMEdge& ced = child_einfo.data;
    bool is_better = false;

    if (use_volume)
    {
      double child_vol = ced.volume(graph.node(child_einfo.source), graph.node(child_einfo.target));
      double current_vol = ced.subtree_volume + child_vol;
      if (current_vol > best_score)
      {
        is_better  = true;
        best_score = current_vol;
      }
    }
    else
    {
      double z         = ced.subtree_max_endZ;
      double child_len = ced.length(graph.node(child_einfo.source), graph.node(child_einfo.target));
      double secondary = ced.subtree_length + child_len;

      if (z > best_score + Z_EPS ||
          (std::abs(z - best_score) <= Z_EPS && secondary > best_secondary))
      {
        is_better      = true;
        best_score     = z;
        best_secondary = secondary;
      }
    }

    if (is_better)
      main_child = child_eid;
  }

  for (int child_eid : child_eids)
  {
    if (child_eid == main_child)
    {
      assign_subtree_ids(child_eid, current_axis_id, current_branch_order, next_axis_id, use_volume);
    }
    else
    {
      int new_id = next_axis_id++;
      assign_subtree_ids(child_eid, new_id, current_branch_order + 1, next_axis_id, use_volume);
    }
  }
}

}

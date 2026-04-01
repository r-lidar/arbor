#include <vector>
#include <algorithm>
#include <map>
#include <limits>
#include <cmath>
#include <unordered_set>
#include <queue>

#include "arbor.h"
#include "QSMbuilder.h"
#include "Grid3D.h"

namespace arbor::qsm {

PointCloud QSMbuilder::clean_tree_butt(const PointCloud& tree)
{
  ServiceLocator::logger()("Checking multiple entry points");

  size_t n = tree.size();
  if (n == 0) return tree;

  double min_z = std::numeric_limits<double>::max();
  for (size_t i = 0; i < n; ++i) {
    if (tree.get_z(i) < min_z) min_z = tree.get_z(i);
  }

  std::vector<bool> is_bottom(n, false);
  std::vector<size_t> bottom_indices;

  for (size_t i = 0; i < n; ++i) {
    if (tree.get_z(i) < (min_z + 0.25)) {
      is_bottom[i] = true;
      bottom_indices.push_back(i);
    }
  }
  PointCloud bottom = tree.subset(is_bottom);

  if (bottom.size() == 0) return tree;

  Grid3D grid(bottom, 0.05);
  std::vector<int> cluster_ids = grid.connected_components(26);

  std::map<int, size_t> counts;
  for (int id : cluster_ids) counts[id]++;

  if (counts.size() <= 1)
    return tree;

  int main_cluster_id = -1;
  size_t max_points = 0;
  for (auto const& [id, count] : counts) {
    if (count > max_points) {
      max_points = count;
      main_cluster_id = id;
    }
  }

  std::vector<bool> to_keep(n, true);
  for (size_t i = 0; i < cluster_ids.size(); ++i)
  {
    if (cluster_ids[i] != main_cluster_id) {
      to_keep[bottom_indices[i]] = false;
    }
  }

  return tree.subset(to_keep);
}

void QSMbuilder::detect_weird_butt(double thresh, int window)
{
  ServiceLocator::logger()("Checking weird butt");

  // Collect main axis edge IDs
  std::vector<int> main_axis_eids;
  for (const auto& [eid, einfo] : graph.edges())
    if (einfo.data.axis_ID == 1) main_axis_eids.push_back(eid);

  if (main_axis_eids.size() < 5) return;

  // Sort root→tip (descending subtree_length)
  std::sort(main_axis_eids.begin(), main_axis_eids.end(), [this](int a, int b) {
    return graph.edge_data(a).subtree_length > graph.edge_data(b).subtree_length;
  });

  size_t i = 0;
  while (i < main_axis_eids.size())
  {
    bool sequence_valid = true;
    for (int w = 0; w < window; ++w)
    {
      size_t idx = i + w;
      if (idx >= main_axis_eids.size())
      {
        sequence_valid = false;
        break;
      }
      int eid = main_axis_eids[idx];
      const auto& einfo = graph.edge(eid);
      double ang = graph.edge_data(eid).angle(
        graph.node(einfo.source),
        graph.node(einfo.target));
      if (ang >= thresh)
      {
        sequence_valid = false;
        break;
      }
    }
    if (sequence_valid) break;
    i++;
  }

  if (i > 0)
  {
    int new_root_eid = (i < main_axis_eids.size()) ? main_axis_eids[i] : -1;

    // Remove weird edges
    for (size_t j = 0; j < i; ++j)
      graph.remove_edge(main_axis_eids[j]);

    remove_disconnected_branches();

    // Set new root's parent_ID to 0
    if (new_root_eid >= 0 && graph.has_edge(new_root_eid))
      graph.edge_data(new_root_eid).parent_ID = 0;
  }
}

void QSMbuilder::remove_disconnected_branches()
{
  // BFS from all main-axis edges to find all reachable edges
  std::unordered_set<int> keep_eids;
  std::queue<int> bfs;

  // Seed BFS with all main axis edges
  for (const auto& [eid, einfo] : graph.edges())
  {
    if (einfo.data.axis_ID == 1)
    {
      keep_eids.insert(eid);
      bfs.push(eid);
    }
  }

  while (!bfs.empty())
  {
    int cur_eid = bfs.front();
    bfs.pop();

    QSM::NodeID tgt = graph.edge(cur_eid).target;
    for (int child_eid : graph.outgoing_edges(tgt))
    {
      if (!keep_eids.count(child_eid))
      {
        keep_eids.insert(child_eid);
        bfs.push(child_eid);
      }
    }
  }

  // Remove edges not in keep set
  std::vector<int> to_remove;
  for (const auto& [eid, einfo] : graph.edges())
    if (!keep_eids.count(eid)) to_remove.push_back(eid);

  for (int eid : to_remove) graph.remove_edge(eid);

  // Remove orphaned nodes (no incident edges)
  std::unordered_set<int> used_nodes;
  for (const auto& [eid, einfo] : graph.edges())
  {
    used_nodes.insert(einfo.source);
    used_nodes.insert(einfo.target);
  }

  std::vector<int> orphan_nodes;
  for (const auto& [nid, ndata] : graph.nodes())
    if (!used_nodes.count(nid)) orphan_nodes.push_back(nid);

  for (int nid : orphan_nodes) graph.remove_node(nid);

  compute_topology();
}

}

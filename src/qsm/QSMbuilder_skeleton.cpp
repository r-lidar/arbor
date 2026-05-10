/**
 * @file QSMbuilder_skeleton.cpp
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

#include <unordered_map>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <algorithm>
#include <sstream>
#include <queue>

#include "arbor.h"
#include "QSMbuilder.h"
#include "fitting.h"
#include "nanoflann.h"

namespace arbor::qsm {

struct ClusterCenter
{
  double x, y, z;
  int iter, id;
  bool done = false;
};

struct pair_hash
{
  std::size_t operator()(const std::pair<int,int>& p) const { return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1); }
};

// Adaptor: nanoflann reads directly from our vector
struct CenterCloud
{
  const std::vector<ClusterCenter>& pts;
  explicit CenterCloud(const std::vector<ClusterCenter>& p) : pts(p) {}
  std::size_t kdtree_get_point_count()               const { return pts.size(); }
  double kdtree_get_pt(std::size_t i, int d)         const { return d == 0 ? pts[i].x : d == 1 ? pts[i].y : pts[i].z; }
  template<class BBOX> bool kdtree_get_bbox(BBOX&)   const { return false; }
};


void QSMbuilder::build_skeleton(const PointCloud& pc, const std::vector<std::pair<int,int>>& iter_cluster, double max_d)
{
  ServiceLocator::logger()("Constructing skeleton");

  // Step 1: group points and compute centers
  // ----------------------------------------
  typedef std::pair<int,int> ClusterKey;
  std::unordered_map<ClusterKey, std::vector<int>, pair_hash> cluster_indices;
  for (size_t i = 0, n = pc.size(); i < n; ++i)
  {
    cluster_indices[iter_cluster[i]].push_back((int)i);
  }

  std::vector<ClusterCenter> centers;
  centers.reserve(cluster_indices.size());
  int id = 1;

  for (auto& [key, indices_binding] : cluster_indices)
  {
    // Re-assign to a standard reference for clang
    auto& indices = indices_binding;

    ClusterCenter c;
    c.iter = key.first;
    c.id   = id++;

    auto compute_mean = [&]()
    {
      double sx = 0, sy = 0, sz = 0;
      for (int idx : indices) { sx += pc.get_x(idx); sy += pc.get_y(idx); sz += pc.get_z(idx); }
      c.x = sx / indices.size();
      c.y = sy / indices.size();
      c.z = sz / indices.size();
    };

    if (indices.size() >= 100)
    {
      utils::fitting::FittingCircloid rc;
      for (int idx : indices)
      {
        rc.add_point(pc.get_x(idx), pc.get_y(idx), pc.get_z(idx));
      }

      utils::fitting::FittingResult ans = rc.fit(0.03);
      if (ans.is_valid(50, 30, 120.0))
      {
        auto& ctr = ans.center; c.x = ctr.x; c.y = ctr.y; c.z = ctr.z;
      }
      else
        compute_mean();
    }
    else
    {
      compute_mean();
    }

    centers.push_back(c);
  }

  if (centers.empty()) return;

  //Step 2: build a static nanoflann KD-tree over all centers
  // --------------------------------------------------------
  using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, CenterCloud>, CenterCloud, 3>;

  CenterCloud cloud(centers);
  KDTree kdtree(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  kdtree.buildIndex();

  const double max_d2 = max_d * max_d;

  // Step 3: min-heap for fast "undone center with smallest iter" lookup
  // (used only in the fallback branch)
  // ---------------------------------

  auto heap_cmp = [](ClusterCenter* a, ClusterCenter* b){ return a->iter > b->iter; };
  std::priority_queue<ClusterCenter*,  std::vector<ClusterCenter*>, decltype(heap_cmp)> minHeap(heap_cmp);

  // Step 4: find initial root (min Z)
  // ---------------------------------
  ClusterCenter* root = &*std::min_element(centers.begin(), centers.end(), [](const ClusterCenter& a, const ClusterCenter& b){ return a.z < b.z; });

  root->done = true;

  for (auto& c : centers)
  {
    if (&c != root) minHeap.push(&c);
  }

  std::unordered_map<int, NodeID> center_to_node;
  center_to_node[root->id] = graph.add_node({root->x, root->y, root->z});

  int remaining = (int)centers.size() - 1;
  int cyl_ID    = 1;

  nanoflann::SearchParameters search_params;
  search_params.sorted = false;

  // Step 5: greedy chain growing loop
  // ---------------------------------
  while (remaining > 0)
  {
    // --- Hot path: radius search around root, filter by iter and done
    double query[3] = { root->x, root->y, root->z };
    std::vector<nanoflann::ResultItem<uint32_t, double>> hits;
    kdtree.radiusSearch(query, max_d2, hits, search_params);

    ClusterCenter* newRoot = nullptr;
    double bestD2 = std::numeric_limits<double>::max();

    for (auto& [idx, d2] : hits)
    {
      ClusterCenter* c = &centers[idx];
      if (c->done || c->iter <= root->iter) continue;
      if (d2 < bestD2) { bestD2 = d2; newRoot = c; }
    }

    if (newRoot)
    {
      // Connect root -> newRoot
      newRoot->done = true;
      remaining--;

      if (!center_to_node.count(root->id))
        center_to_node[root->id] = graph.add_node({root->x, root->y, root->z});
      if (!center_to_node.count(newRoot->id))
        center_to_node[newRoot->id] = graph.add_node({newRoot->x, newRoot->y, newRoot->z});

      QSMEdge ed; ed.cyl_ID = cyl_ID++;
      graph.add_edge(center_to_node[root->id], center_to_node[newRoot->id], ed);

      root = newRoot;
    }
    else
    {
      // Fallback: chain is stuck, start a new branch

      // Drain stale (already-done) entries from the heap - O(log n) amortised
      while (!minHeap.empty() && minHeap.top()->done) minHeap.pop();
      if (minHeap.empty()) break;

      ClusterCenter* orphan = minHeap.top();
      minHeap.pop();

      // Find nearest already-done center to orphan (O(n_done), infrequent)
      // This scan fires only O(branch_count) times, not O(n²)
      ClusterCenter* nearestDone = nullptr;
      double bestDist = std::numeric_limits<double>::max();
      for (auto& c : centers)
      {
        if (!c.done) continue;
        double dx = c.x - orphan->x, dy = c.y - orphan->y, dz = c.z - orphan->z;
        double d  = dx*dx + dy*dy + dz*dz;
        if (d < bestDist) { bestDist = d; nearestDone = &c; }
      }

      if (!nearestDone) break;

      orphan->done = true;
      remaining--;

      if (!center_to_node.count(nearestDone->id))
        center_to_node[nearestDone->id] = graph.add_node({nearestDone->x, nearestDone->y, nearestDone->z});
      if (!center_to_node.count(orphan->id))
        center_to_node[orphan->id] = graph.add_node({orphan->x, orphan->y, orphan->z});

      QSMEdge ed; ed.cyl_ID = cyl_ID++;
      graph.add_edge(center_to_node[nearestDone->id], center_to_node[orphan->id], ed);

      root = orphan;
    }
  }
}

void QSMbuilder::fix_multiple_root()
{
  NodeID old_root_id = -1;

  // Find the single existing root
  for (const auto& [nid, _] : graph.nodes())
  {
    if (graph.incoming_edges(nid).empty())
    {
      old_root_id = nid;
      break;
    }
  }

  if (old_root_id == -1) return; // No nodes exist yet

  // Create the new root 1mm below
  const QSMNode& old_root = graph.node(old_root_id);
  QSMNode new_root_node = old_root;
  new_root_node.z -= 0.01;

  NodeID new_root_id = graph.add_node(new_root_node);

  // Connect them
  QSMEdge new_edge_data;
  new_edge_data.cyl_ID = static_cast<int>(graph.edge_count()) + 1;

  graph.add_edge(new_root_id, old_root_id, new_edge_data);

  compute_topology();
}

static void collect_subtree(QSM& graph, NodeID  node_id, std::vector<EdgeID>& edges_out, std::vector<NodeID>& nodes_out)
{
  for (EdgeID eid : graph.outgoing_edges(node_id))
  {
    NodeID child = graph.edge(eid).target;
    edges_out.push_back(eid);
    nodes_out.push_back(child);
    collect_subtree(graph, child, edges_out, nodes_out);
  }
}

void QSMbuilder::prune_spurious_branches()
{
  std::unordered_map<int, std::vector<EdgeID>> axis_edges;

  for (const auto& [eid, einfo] : graph.edges())
  {
    if (einfo.data.branch_order == 2)
      axis_edges[einfo.data.axis_ID].push_back(eid);
  }

  std::vector<EdgeID> edges_to_remove;
  std::vector<NodeID> nodes_to_remove;

  for (const auto& [axis_id, eids] : axis_edges)
  {
    if (eids.size() != 1) continue;

    EdgeID eid = eids[0];
    NodeID subtree_root = graph.edge(eid).target;

    edges_to_remove.push_back(eid);
    nodes_to_remove.push_back(subtree_root);
    collect_subtree(graph, subtree_root, edges_to_remove, nodes_to_remove);
  }

  for (EdgeID eid : edges_to_remove) graph.remove_edge(eid);
  for (NodeID nid : nodes_to_remove) graph.remove_node(nid);
}

}

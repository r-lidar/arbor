#include <unordered_map>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <algorithm>
#include <sstream>

#include "arbor.h"
#include "QSMbuilder.h"
#include "ransac.h"

namespace arbor::qsm {

void QSMbuilder::build_skeleton(const PointCloud& pc, const std::vector<std::pair<int, int>>& iter_cluster, double max_d)
{
  struct ClusterCenter {
    double x, y, z;
    int iter;
    int id;
    bool done = false;
  };

  struct pair_hash {
    std::size_t operator()(const std::pair<int, int>& p) const {
      return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
  };

  logger("Constructing skeleton");

  size_t n = pc.size();

  // Step 1: compute mean for each (iter, cluster)
  typedef std::pair<int, int> ClusterKey;
  std::unordered_map<ClusterKey, std::vector<int>, pair_hash> cluster_indices;

  for (size_t i = 0; i < n; ++i) {
    cluster_indices[iter_cluster[i]].push_back(i);
  }

  // build centers using RANSAC
  std::vector<ClusterCenter> centers;
  int id = 1;

  for (auto& entry : cluster_indices)
  {
    auto key = entry.first;
    auto& indices = entry.second;

    ClusterCenter c;
    c.iter = key.first;
    c.id = id++;

    if (indices.size() >= 100)
    {
      RansacCircle rc(100, 0.02);

      for (int idx : indices)
      {
        double x = pc.get_x(idx);
        double y = pc.get_y(idx);
        double z = pc.get_z(idx);
        rc.add_point(x, y, z);
      }

      rc.find_circle();

      if (rc.is_valid(0.5, 0.3, 120.0))
      {
        auto center = rc.get_center();
        c.x = center[0];
        c.y = center[1];
        c.z = center[2];
      }
      else
      {
        double sumx = 0, sumy = 0, sumz = 0;
        for (int idx : indices)
        {
          sumx += pc.get_x(idx);
          sumy += pc.get_y(idx);
          sumz += pc.get_z(idx);
        }
        c.x = sumx / indices.size();
        c.y = sumy / indices.size();
        c.z = sumz / indices.size();
      }
    }
    else
    {
      double sumx = 0, sumy = 0, sumz = 0;
      for (int idx : indices)
      {
        sumx += pc.get_x(idx);
        sumy += pc.get_y(idx);
        sumz += pc.get_z(idx);
      }
      c.x = sumx / indices.size();
      c.y = sumy / indices.size();
      c.z = sumz / indices.size();
    }

    centers.push_back(c);
  }

  // Prepare for neighbor searching.
  std::vector<ClusterCenter*> searchSpace;
  for (auto& c : centers)
    searchSpace.push_back(&c);

  if (searchSpace.empty()) return;

  // Find initial root: cluster center with minimum Z value.
  ClusterCenter* root = *std::min_element(
    searchSpace.begin(),
    searchSpace.end(),
    [](ClusterCenter* a, ClusterCenter* b) { return a->z < b->z; }
  );

  root->done = true;
  searchSpace.erase(std::remove(searchSpace.begin(), searchSpace.end(), root), searchSpace.end());

  // Map ClusterCenter id -> graph NodeID
  std::unordered_map<int, QSMGraph::NodeID> center_to_node;
  center_to_node[root->id] = graph.add_node({root->x, root->y, root->z});

  const double max_d2 = max_d * max_d;

  int cyl_ID = 1;

  while (!searchSpace.empty())
  {
    ClusterCenter* start   = root;
    ClusterCenter* newRoot = nullptr;
    double bestD2 = std::numeric_limits<double>::max();

    for (auto* c : searchSpace)
    {
      if (c->iter <= root->iter) continue;

      double dx = c->x - root->x;
      double dy = c->y - root->y;
      double dz = c->z - root->z;
      double d2 = dx*dx + dy*dy + dz*dz;

      if (d2 < bestD2 && d2 <= max_d2)
      {
        bestD2 = d2;
        newRoot = c;
      }
    }

    if (newRoot)
    {
      newRoot->done = true;
      searchSpace.erase(std::remove(searchSpace.begin(), searchSpace.end(), newRoot), searchSpace.end());

      if (!center_to_node.count(start->id))
        center_to_node[start->id] = graph.add_node({start->x, start->y, start->z});
      if (!center_to_node.count(newRoot->id))
        center_to_node[newRoot->id] = graph.add_node({newRoot->x, newRoot->y, newRoot->z});

      QSMEdge ed;
      ed.cyl_ID = cyl_ID++;
      graph.add_edge(center_to_node[start->id], center_to_node[newRoot->id], ed);

      root = newRoot;
    }
    else
    {
      auto minIt = std::min_element(
        searchSpace.begin(),
        searchSpace.end(),
        [](ClusterCenter* a, ClusterCenter* b) { return a->iter < b->iter; }
      );
      root = *minIt;

      ClusterCenter* nearestDone = nullptr;
      double bestDist = std::numeric_limits<double>::max();
      for (auto& c : centers)
      {
        if (!c.done) continue;

        double dx = c.x - root->x;
        double dy = c.y - root->y;
        double dz = c.z - root->z;
        double d = dx*dx + dy*dy + dz*dz;

        if (d < bestDist)
        {
          bestDist = d;
          nearestDone = &c;
        }
      }

      if (!nearestDone) break;

      root->done = true;
      searchSpace.erase(minIt);

      if (!center_to_node.count(nearestDone->id))
        center_to_node[nearestDone->id] = graph.add_node({nearestDone->x, nearestDone->y, nearestDone->z});
      if (!center_to_node.count(root->id))
        center_to_node[root->id] = graph.add_node({root->x, root->y, root->z});

      QSMEdge ed;
      ed.cyl_ID = cyl_ID++;
      graph.add_edge(center_to_node[nearestDone->id], center_to_node[root->id], ed);
    }
  }
}

void QSMbuilder::fix_multiple_root()
{
  QSMGraph::NodeID old_root_id = -1;

  // Find the single existing root
  for (const auto& [nid, _] : graph.nodes())
  {
    if (graph.incoming_edges(nid).empty())
    {
      old_root_id = nid;
      break;
    }
  }

  printf("old_root_id %d\n", old_root_id);

  if (old_root_id == -1) return; // No nodes exist yet

  // Create the new root 1mm below
  const QSMNode& old_root = graph.node(old_root_id);
  QSMNode new_root_node = old_root;
  new_root_node.z -= 0.01;

  QSMGraph::NodeID new_root_id = graph.add_node(new_root_node);

  // Connect them
  QSMEdge new_edge_data;
  new_edge_data.cyl_ID = static_cast<int>(graph.edge_count()) + 1;

  graph.add_edge(new_root_id, old_root_id, new_edge_data);

  compute_topology();
}

}

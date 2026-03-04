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

PointCloud QSMbuilder::clean_tree_butt(const PointCloud& tree, const Logger& logger)
{
  logger("Checking multiple entry points");

  size_t n = tree.size();
  if (n == 0) return tree;

  // 1. Find Min Z
  double min_z = std::numeric_limits<double>::max();
  for (size_t i = 0; i < n; ++i) {
    if (tree.get_z(i) < min_z) min_z = tree.get_z(i);
  }

  // 2. Identify "bottom" points and keep track of original indices
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

  // 3. Connected Components
  // Parameters from R: resolution = 0.05, connectivity = 26
  Grid3D grid(bottom, 0.05);
  std::vector<int> cluster_ids = grid.connected_components(26);

  // 4. Analyze clusters
  std::map<int, size_t> counts;
  for (int id : cluster_ids) {
    counts[id]++;
  }

  // If more than one cluster is detected
  if (counts.size() <= 1)
    return tree;

  // Find the ID of the largest cluster (the main trunk)
  int main_cluster_id = -1;
  size_t max_points = 0;
  for (auto const& [id, count] : counts) {
    if (count > max_points) {
      max_points = count;
      main_cluster_id = id;
    }
  }

  // 5. Prepare final subset (remove points from smaller clusters)
  std::vector<bool> to_keep(n, true);
  for (size_t i = 0; i < cluster_ids.size(); ++i)
  {
    // If this point in 'bottom' is NOT part of the main cluster
    if (cluster_ids[i] != main_cluster_id) {
      // Mark the corresponding original index for removal
      size_t original_idx = bottom_indices[i];
      to_keep[original_idx] = false;
    }
  }

  return tree.subset(to_keep);
}

void QSMbuilder::detect_weird_butt(double thresh, int window)
{
  logger("Checking weird butt");

  // Get the main axis (trunk)
  // Store IDs instead of pointers to be safe from map reallocations
  std::vector<int> main_axis_ids;
  for (auto& [id, cyl] : qsm)
  {
    if (cyl.axis_ID == 1) {
      main_axis_ids.push_back(id);
    }
  }

  if (main_axis_ids.empty()) return;

  // Sort IDs based on the subtree_length of the cylinders they represent
  // to get the root at index 0 and then subsequent cylinders
  std::sort(main_axis_ids.begin(), main_axis_ids.end(), [this](int a, int b) { return qsm.cylinders_[a].subtree_length > qsm.cylinders_[b].subtree_length; });

  size_t i = 0;
  while (i < main_axis_ids.size())
  {
    bool sequence_valid = true;
    for (int w = 0; w < window; ++w)
    {
      size_t idx = i + w;
      if (idx >= main_axis_ids.size() || qsm.cylinders_[main_axis_ids[idx]].angle() >= thresh)
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
    // Determine the ID of the new root BEFORE erasing anything
    // The new root is the first valid cylinder in the sequence (index i)
    int new_root_id = (i < main_axis_ids.size()) ? main_axis_ids[i] : -1;

    // Remove the weird cylinders (0 to i-1)
    for (size_t j = 0; j < i; ++j) {
      qsm.cylinders_.erase(main_axis_ids[j]);
    }

    remove_disconnected_branches();

    // Update the new root's parent in the NEW map
    if (new_root_id != -1 && qsm.cylinders_.count(new_root_id))
    {
      qsm.cylinders_[new_root_id].parent_ID = 0;
    }
  }
}

void QSMbuilder::remove_disconnected_branches()
{
  std::unordered_set<int> keep_set;
  std::queue<int> traversal_queue;

  // Start with all cylinders currently on the main axis
  for (auto const& [id, cyl] : qsm)
  {
    if (cyl.axis_ID == 1)
    {
      keep_set.insert(id);
      traversal_queue.push(id);
    }
  }

  // BFS to find all descendants (connected cylinders)
  while (!traversal_queue.empty())
  {
    int current_id = traversal_queue.front();
    traversal_queue.pop();

    if (qsm.children_map_.count(current_id))
    {
      for (int child_id : qsm.children_map_.at(current_id))
      {
        if (keep_set.find(child_id) == keep_set.end())
        {
          keep_set.insert(child_id);
          traversal_queue.push(child_id);
        }
      }
    }
  }

  // Rebuild the cylinder map containing only connected components
  std::unordered_map<int, QSMcylinder> new_cylinders;
  for (int id : keep_set)
  {
    new_cylinders[id] = qsm.cylinders_[id];
  }

  qsm.cylinders_ = std::move(new_cylinders);

  // Refresh topology
  compute_topology();
}

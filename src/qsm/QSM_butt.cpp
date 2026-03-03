#include <vector>
#include <algorithm>
#include <map>
#include <limits>
#include <cmath>
#include <unordered_set>
#include <queue>

#include "arbor.h"
#include "Grid3D.h"

PointCloud QSM::clean_tree_butt(const PointCloud& tree)
{
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

void QSM::detect_weird_butt(double thresh, int window)
{
  std::cout << "Validating butt architecture..." << std::endl;

  // Get the main axis cylinders
  std::vector<QSMcylinder*> main_axis_cyls;
  for (auto& [id, cyl] : cylinders_)
  {
    if (cyl.axis_ID == 1)
    {
      main_axis_cyls.push_back(&cyl);
    }
  }

  if (main_axis_cyls.empty()) return;

  // Sort from root to tip (using subtree_length as a proxy for height/order).
  std::sort(main_axis_cyls.begin(), main_axis_cyls.end(), [](QSMcylinder* a, QSMcylinder* b) { return a->subtree_length > b->subtree_length; });

  // Loop from root. Look a n consecutive cylinders. If they all have an angle < threshold
  // then it is valid. Otherwise we removes cylinder with an angle too steep. They are likely
  // to be artifacts of bad segmentation
  size_t i = 0;
  bool fix_needed = false;

  while (i < main_axis_cyls.size())
  {
    bool sequence_valid = true;

    for (int w = 0; w < window; ++w)
    {
      size_t idx = i + w;
      if (idx >= main_axis_cyls.size())
      {
        sequence_valid = false;
        break;
      }
      if (main_axis_cyls[idx]->angle() >= thresh)
      {
        sequence_valid = false;
        break;
      }
    }

    if (sequence_valid) break;
    i++;
  }

  // 3. If i > 0, we found "weird" cylinders at the start
  if (i > 0)
  {
    std::cerr << "[WARN] Detection of weird tree butt. Automatic fix triggered." << std::endl;

    // Remove the weird cylinders from the map
    for (size_t j = 0; j <= i; ++j)
    {
      cylinders_.erase(main_axis_cyls[j]->cyl_ID);
    }

    // Remove branches that are now floating
    remove_disconnected_branches();

    // Find the new root (the first cylinder of axis 1 remaining)
    main_axis_cyls[0]->parent_ID = 0;
  }
}

void QSM::remove_disconnected_branches()
{
  std::unordered_set<int> keep_set;
  std::queue<int> traversal_queue;

  // Step 1: Start with all cylinders currently on the main axis
  for (auto const& [id, cyl] : cylinders_) {
    if (cyl.axis_ID == 1) {
      keep_set.insert(id);
      traversal_queue.push(id);
    }
  }

  // Step 2: BFS to find all descendants (connected cylinders)
  while (!traversal_queue.empty()) {
    int current_id = traversal_queue.front();
    traversal_queue.pop();

    if (children_map_.count(current_id)) {
      for (int child_id : children_map_.at(current_id)) {
        if (keep_set.find(child_id) == keep_set.end()) {
          keep_set.insert(child_id);
          traversal_queue.push(child_id);
        }
      }
    }
  }

  // Step 3: Rebuild the cylinder map containing only connected components
  std::unordered_map<int, QSMcylinder> new_cylinders;
  for (int id : keep_set) {
    new_cylinders[id] = cylinders_[id];
  }

  cylinders_ = std::move(new_cylinders);

  // Refresh topology
  compute_topology();
}

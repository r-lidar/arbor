#include "QSM.h"
#include <algorithm>
#include <cmath>
#include <limits>

// small epsilon for z comparisons
static constexpr double Z_EPS = 1e-9;

void QSM::build_from_vectors(const std::vector<int>& cyl_ID,
                             const std::vector<int>& parent_ID,
                             const std::vector<double>& length,
                             const std::vector<double>& startZ,
                             const std::vector<double>& endZ)
{
  children_map_.clear();
  length_map_.clear();
  startz_map_.clear();
  endz_map_.clear();

  size_t n = cyl_ID.size();
  for (size_t i = 0; i < n; ++i)
  {
    int cid = cyl_ID[i];
    int pid = parent_ID[i];

    // note: allow multiple children, parent may be 0 or -1 if root's parent
    children_map_[pid].push_back(cid);

    // ensure node exists in children map (so leaf nodes have entry only if referenced)
    if (children_map_.find(cid) == children_map_.end())
      children_map_.emplace(cid, std::vector<int>());

    length_map_[cid] = length[i];
    startz_map_[cid] = startZ[i];
    endz_map_[cid] = endZ[i];
  }
}

double QSM::compute_subtree_length(NodeId node_id)
{
  auto itCached = subtree_lengths_.find(node_id);
  if (itCached != subtree_lengths_.end())
    return itCached->second;

  auto it = children_map_.find(node_id);
  if (it == children_map_.end() || it->second.empty())
  {
    subtree_lengths_[node_id] = 0.0;
    return 0.0;
  }

  double max_length = 0.0;
  for (NodeId child_id : it->second)
  {
    double child_length = compute_subtree_length(child_id);
    // length_map_ contains child own length (length of the cylinder from child)
    double candidate = child_length + (length_map_.count(child_id) ? length_map_.at(child_id) : 0.0);
    if (candidate > max_length)
      max_length = candidate;
  }

  subtree_lengths_[node_id] = max_length;
  return max_length;
}

double QSM::compute_subtree_max_z(NodeId node_id)
{
  auto itCached = subtree_max_z_.find(node_id);
  if (itCached != subtree_max_z_.end())
    return itCached->second;

  // start with this node's own endZ (it may be a tip)
  double my_max_z = std::numeric_limits<double>::lowest();
  if (endz_map_.count(node_id))
    my_max_z = endz_map_.at(node_id);
  else
    my_max_z = std::numeric_limits<double>::lowest();

  auto it = children_map_.find(node_id);
  if (it != children_map_.end())
  {
    for (NodeId child_id : it->second)
    {
      double child_max_z = compute_subtree_max_z(child_id);
      if (child_max_z > my_max_z)
        my_max_z = child_max_z;
    }
  }

  subtree_max_z_[node_id] = my_max_z;
  return my_max_z;
}

void QSM::assign_subtree_ids(NodeId node_id, int current_subtree_id, int current_branching_order, int &next_subtree_id)
{
  subtree_ids_[node_id] = current_subtree_id;
  branching_orders_[node_id] = current_branching_order;

  auto it = children_map_.find(node_id);
  if (it == children_map_.end() || it->second.empty())
    return;

  // choose main child by subtree_max_z (highest tip). tie-breaker: subtree length + length(child)
  NodeId main_child = -1;
  double best_max_z = -std::numeric_limits<double>::infinity();
  double best_secondary = -std::numeric_limits<double>::infinity();

  for (NodeId child_id : it->second)
  {
    double child_max_z = subtree_max_z_.count(child_id) ? subtree_max_z_.at(child_id)
      : (endz_map_.count(child_id) ? endz_map_.at(child_id) : -std::numeric_limits<double>::infinity());
    double secondary = 0.0;
    if (subtree_lengths_.count(child_id))
      secondary = subtree_lengths_.at(child_id) + (length_map_.count(child_id) ? length_map_.at(child_id) : 0.0);

    if ( (child_max_z > best_max_z + Z_EPS) ||
         (std::fabs(child_max_z - best_max_z) <= Z_EPS && secondary > best_secondary) )
    {
      best_max_z = child_max_z;
      best_secondary = secondary;
      main_child = child_id;
    }
  }

  for (NodeId child_id : it->second)
  {
    if (child_id == main_child)
    {
      assign_subtree_ids(child_id, current_subtree_id, current_branching_order, next_subtree_id);
    }
    else
    {
      int new_subtree_id = next_subtree_id++;
      assign_subtree_ids(child_id, new_subtree_id, current_branching_order + 1, next_subtree_id);
    }
  }
}

void QSM::compute_architecture(NodeId root_id)
{
  // clear previous results
  subtree_lengths_.clear();
  subtree_max_z_.clear();
  subtree_ids_.clear();
  branching_orders_.clear();

  // compute subtree lengths starting from root (post-order recursion)
  compute_subtree_length(root_id);

  // compute subtree max z
  compute_subtree_max_z(root_id);

  // assign subtree ids, starting with axis id 1 and branching order 1
  int next_subtree_id = 2;
  assign_subtree_ids(root_id, 1, 1, next_subtree_id);
}

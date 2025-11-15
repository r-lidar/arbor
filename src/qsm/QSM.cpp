#include "QSM.h"

#include <algorithm>
#include <cmath>
#include <limits>

void QSM::build_from_cylinders(const std::vector<QSMcylinder>& input)
{
  cylinders_.clear();
  children_map_.clear();

  // Insert cylinders and build children links
  for (const auto& c : input)
  {
    cylinders_[c.cyl_ID] = c;

    // Ensure parent exists in the child map
    children_map_[c.parent_ID].push_back(c.cyl_ID);

    // Ensure the child also has a children entry even if empty
    if (!children_map_.count(c.cyl_ID))
      children_map_[c.cyl_ID] = {};
  }
}

double QSM::compute_subtree_length(int node_id)
{
  auto& node = cylinders_[node_id];

  // Cache hit?
  if (node.subtree_length >= 0)
    return node.subtree_length;

  const auto& kids = children_map_[node_id];
  if (kids.empty())
  {
    node.subtree_length = 0.0;
    return 0.0;
  }

  double max_len = 0.0;
  for (int child_id : kids)
  {
    auto& child = cylinders_[child_id];
    double candidate = compute_subtree_length(child_id) + child.length();
    max_len = std::max(max_len, candidate);
  }

  node.subtree_length = max_len;
  return max_len;
}

double QSM::compute_subtree_max_z(int node_id)
{
  auto& node = cylinders_[node_id];

  // Cache hit?
  if (node.subtree_max_endZ > SUBTREE_MAXZ_UNSET)
    return node.subtree_max_endZ;

  double maxz = node.endZ;

  for (int child_id : children_map_[node_id])
  {
    double child_maxz = compute_subtree_max_z(child_id);
    maxz = std::max(maxz, child_maxz);
  }

  node.subtree_max_endZ = maxz;
  return maxz;
}

void QSM::assign_subtree_ids(int node_id, int current_axis_id, int current_branch_order, int &next_axis_id)
{
  auto& node = cylinders_[node_id];
  node.axis_ID = current_axis_id;
  node.branch_order = current_branch_order;

  const auto& kids = children_map_[node_id];
  if (kids.empty()) return;

  // Select main child: highest endZ, then longest subtree
  int main_child = -1;
  double bestZ = -1e300;
  double bestSecondary = -1e300;

  for (int child_id : kids)
  {
    auto& child = cylinders_[child_id];

    double z = child.subtree_max_endZ;
    double secondary = child.subtree_length + child.length();

    if (z > bestZ + Z_EPS || (std::abs(z - bestZ) <= Z_EPS && secondary > bestSecondary))
    {
      main_child = child_id;
      bestZ = z;
      bestSecondary = secondary;
    }
  }

  for (int child_id : kids)
  {
    if (child_id == main_child)
    {
      assign_subtree_ids(child_id, current_axis_id, current_branch_order, next_axis_id);
    }
    else
    {
      int new_id = next_axis_id++;
      assign_subtree_ids(child_id, new_id, current_branch_order + 1, next_axis_id);
    }
  }
}

void QSM::compute_architecture(int root_id)
{
  // initialize caches inside cylinders
  for (auto& kv : cylinders_)
  {
    kv.second.subtree_length    = SUBTREE_LENGTH_UNSET;
    kv.second.subtree_max_endZ  = SUBTREE_MAXZ_UNSET;
    kv.second.axis_ID           = 0;
    kv.second.branch_order      = 0;
  }

  compute_subtree_length(root_id);
  compute_subtree_max_z(root_id);

  int next_axis_id = 2;
  assign_subtree_ids(root_id, 1, 1, next_axis_id);
}

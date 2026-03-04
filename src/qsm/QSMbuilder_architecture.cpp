#include "QSMbuilder.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

void QSMbuilder::compute_topology()
{
  logger("Connecting nodes");

  struct CoordKey {
    int x, y, z;
    bool operator==(const CoordKey& other) const noexcept {
      return x == other.x && y == other.y && z == other.z;
    }
  };

  // 2. Define local hasher
  struct CoordKeyHash {
    std::size_t operator()(const CoordKey& k) const noexcept {
      std::size_t h1 = std::hash<int>{}(k.x);
      std::size_t h2 = std::hash<int>{}(k.y);
      std::size_t h3 = std::hash<int>{}(k.z);
      return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
  };

  auto make_coord_key = [](double x, double y, double z, int digits = 6) {
    double factor = std::pow(10.0, digits);
    return CoordKey{
      static_cast<int>(std::llround(x * factor)),
      static_cast<int>(std::llround(y * factor)),
      static_cast<int>(std::llround(z * factor))
    };
  };

  qsm.children_map_.clear();

  // Build lookup: end_key -> cyl_ID
  std::unordered_map<CoordKey, int, CoordKeyHash> end_lookup;
  end_lookup.reserve(qsm.cylinders_.size());

  for (auto& [cid, cyl] : qsm)
  {
    CoordKey key = make_coord_key(cyl.endX, cyl.endY, cyl.endZ);
    end_lookup.emplace(key, cid);
  }

  // Now assign parent_ID by matching start_key to end_key
  for (auto& [cid, cyl] : qsm)
  {
    CoordKey start_key = make_coord_key(cyl.startX, cyl.startY, cyl.startZ);

    auto it = end_lookup.find(start_key);
    if (it != end_lookup.end())
    {
      cyl.parent_ID = it->second;
      qsm.children_map_[it->second].push_back(cid);
    }
    else
    {
      cyl.parent_ID = 0;  // no parent -> root
    }
  }
}

void QSMbuilder::compute_architecture(int root_id, bool use_volume)
{
  logger("Computing architecture");

  // Initialize caches inside cylinders
  for (auto& kv : qsm)
  {
    kv.second.subtree_length    = SUBTREE_LENGTH_UNSET;
    kv.second.subtree_max_endZ  = SUBTREE_MAXZ_UNSET;
    kv.second.subtree_volume    = SUBTREE_VOLUME_UNSET;
    kv.second.axis_ID           = 0;
    kv.second.branch_order      = 0;
  }

  compute_subtree_length(root_id);

  if (use_volume)
    compute_subtree_volume(root_id);
  else
    compute_subtree_max_z(root_id);

  int next_axis_id = 2;
  assign_subtree_ids(root_id, 1, 1, next_axis_id, use_volume);
}

double QSMbuilder::compute_subtree_length(int node_id)
{
  auto& node = qsm.cylinders_[node_id];

  // Cache hit?
  if (node.subtree_length >= 0)
    return node.subtree_length;

  const auto& kids = qsm.children_map_[node_id];
  if (kids.empty())
  {
    node.subtree_length = 0.0;
    return 0.0;
  }

  double max_len = 0.0;
  for (int child_id : kids)
  {
    auto& child = qsm.cylinders_[child_id];
    double candidate = compute_subtree_length(child_id) + child.length();
    max_len = std::max(max_len, candidate);
  }

  node.subtree_length = max_len;
  return max_len;
}

double QSMbuilder::compute_subtree_max_z(int node_id)
{
  auto& node = qsm.cylinders_[node_id];

  // Cache hit?
  if (node.subtree_max_endZ > SUBTREE_MAXZ_UNSET)
    return node.subtree_max_endZ;

  double maxz = node.endZ;

  for (int child_id : qsm.children_map_[node_id])
  {
    double child_maxz = compute_subtree_max_z(child_id);
    maxz = std::max(maxz, child_maxz);
  }

  node.subtree_max_endZ = maxz;
  return maxz;
}

double QSMbuilder::compute_subtree_volume(int node_id)
{
  auto& node = qsm.cylinders_[node_id];

  // Cache hit?
  if (node.subtree_volume >= 0)
    return node.subtree_volume;

  const auto& kids = qsm.children_map_[node_id];
  if (kids.empty())
  {
    node.subtree_volume = 0.0;
    return 0.0;
  }

  double max_v = 0.0;
  for (int child_id : kids)
  {
    auto& child = qsm.cylinders_[child_id];
    double candidate = compute_subtree_volume(child_id) + child.volume();
    max_v = std::max(max_v, candidate);
  }

  node.subtree_length = max_v;
  return max_v;
}

void QSMbuilder::assign_subtree_ids(int node_id, int current_axis_id, int current_branch_order, int &next_axis_id, bool use_volume)
{
  auto& node = qsm.cylinders_[node_id];
  node.axis_ID = current_axis_id;
  node.branch_order = current_branch_order;

  const auto& kids = qsm.children_map_[node_id];
  if (kids.empty()) return;

  int main_child = -1;
  double best_score = -1e300;     // General score variable
  double best_secondary = -1e300; // Used only for Z logic tie-breaking

  for (int child_id : kids)
  {
    auto& child = qsm.cylinders_[child_id];
    bool is_better = false;

    if (use_volume)
    {
      // Strategy A: Follow the path with the Highest Volume
      // We look at the child's path volume + the child's own volume
      double current_vol = child.subtree_volume + child.volume();
      if (current_vol > best_score)
      {
        is_better = true;
        best_score = current_vol;
      }
    }
    else
    {
      // Strategy B: Original (Height + Length tie-breaker)
      double z = child.subtree_max_endZ;
      double secondary = child.subtree_length + child.length();

      if (z > best_score + Z_EPS || (std::abs(z - best_score) <= Z_EPS && secondary > best_secondary))
      {
        is_better = true;
        best_score = z;
        best_secondary = secondary;
      }
    }

    if (is_better)
    {
      main_child = child_id;
    }
  }

  for (int child_id : kids)
  {
    if (child_id == main_child)
    {
      // Continue the current axis
      assign_subtree_ids(child_id, current_axis_id, current_branch_order, next_axis_id, use_volume);
    }
    else
    {
      // Start a new axis
      int new_id = next_axis_id++;
      assign_subtree_ids(child_id, new_id, current_branch_order + 1, next_axis_id, use_volume);
    }
  }
}

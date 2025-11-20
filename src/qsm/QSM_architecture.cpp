#include "QSM.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

struct CoordKey
{
  int x, y, z;
  bool operator==(const CoordKey& other) const noexcept { return x == other.x && y == other.y && z == other.z; }
};

struct CoordKeyHash
{
  std::size_t operator()(const CoordKey& k) const noexcept
  {
    std::size_t h1 = std::hash<int>{}(k.x);
    std::size_t h2 = std::hash<int>{}(k.y);
    std::size_t h3 = std::hash<int>{}(k.z);
    return h1 ^ (h2 << 1) ^ (h3 << 2);
  }
};

inline CoordKey make_coord_key(double x, double y, double z, int digits = 6)
{
  double factor = std::pow(10.0, digits);
  return CoordKey{
    static_cast<int>(std::llround(x * factor)),
    static_cast<int>(std::llround(y * factor)),
    static_cast<int>(std::llround(z * factor))
  };
}

void QSM::compute_topology()
{
  children_map_.clear();

  // Build lookup: end_key -> cyl_ID
  std::unordered_map<CoordKey, int, CoordKeyHash> end_lookup;
  end_lookup.reserve(cylinders_.size());

  for (auto& [cid, cyl] : cylinders_)
  {
    CoordKey key = make_coord_key(cyl.endX, cyl.endY, cyl.endZ);
    end_lookup.emplace(key, cid);
  }

  // Now assign parent_ID by matching start_key to end_key
  for (auto& [cid, cyl] : cylinders_)
  {
    CoordKey start_key = make_coord_key(cyl.startX, cyl.startY, cyl.startZ);

    auto it = end_lookup.find(start_key);
    if (it != end_lookup.end())
    {
      cyl.parent_ID = it->second;
      children_map_[it->second].push_back(cid);
    }
    else
    {
      cyl.parent_ID = 0;  // no parent → root
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

// Returns the ordered main axis as pointers (root -> tip)
// axis_ID == 1
std::vector<const QSMcylinder*> QSM::main_axis() const
{
  std::vector<const QSMcylinder*> axis;
  axis.reserve(cylinders_.size());

  for (const auto &kv : cylinders_)
  {
    const QSMcylinder* c = &kv.second;
    if (c->axis_ID == 1)
    {
      axis.push_back(c);
    }
  }

  std::sort(axis.begin(), axis.end(),  [](const QSMcylinder* a, const QSMcylinder* b)
  {
    return a->cyl_ID < b->cyl_ID;
  });

  return axis;
}

void QSM::prolongate(double d, double L)
{
  if (d <= 0.0) return;

  // 1. Get ordered main axis
  std::vector<const QSMcylinder*> axis = main_axis();
  const size_t n = axis.size();

  if (n < 2) return;

  // 2. Compute cyl lengths and cumulative lengths
  std::vector<double> lens(n);
  std::vector<double> cum(n);

  double total = 0.0;
  for (size_t i = 0; i < n; ++i)
  {
    lens[i] = axis[i]->length();
    total += lens[i];
    cum[i] = total;
  }
  if (total == 0.0) return;

  // 3. take only few first cylinders
  // first 10% of axis or 30 cm
  double cutoff = 0.1 * total;
  size_t k = 0;
  while (k < n && cum[k] <= cutoff) k++;
  if (k == 0) k = 1; // minimal 1 cyl
  if (k < n)
  {
    // Ensure sum < 0.3 : add next cylinder
    double s = 0.0;
    for (size_t i = 0; i < k; i++) s += lens[i];
    if (s < 0.3) k++;
  }
  if (k > n) k = n;

  const QSMcylinder* root = axis[0];
  const QSMcylinder* last  = axis[k-1];

  if (root->parent_ID != 0)
    throw std::runtime_error("Invalid QSM, the selected root parent ID is not 0");

  // Orientation estimated
  double dx = last->endX - root->startX;
  double dy = last->endY - root->startY;
  double dz = last->endZ - root->startZ;
  double N  = std::sqrt(dx*dx + dy*dy + dz*dz);
  if (N <= 0.0) return;
  double ox = dx / N;
  double oy = dy / N;
  double oz = dz / N;

  // Adjust prolongation distance by angle
  double d_adj;
  if (std::abs(oz) < 1e-9)
    d_adj = d;
  else
    d_adj = d / oz;

  // Start/end points of prolongation is the root
  double endX = root->startX;
  double endY = root->startY;
  double endZ = root->startZ;

  // Create new subdivided segments
  int nseg = std::max(1, int(std::ceil(d_adj / L)));
  double actual_L = d_adj / nseg;

  double root_subtree = root->subtree_length;
  if (root_subtree == SUBTREE_LENGTH_UNSET)
    throw std::runtime_error("Invalid QSM, the root has not subtree length");

  // Find safe negative cyl_ID
  int next_id = -1;
  int prev_id = 0;

  for (int i = 1; i <= nseg; i++)
  {
    double f1 = double(i - 1) / nseg;
    double f2 = double(i)     / nseg;

    double x1 = endX - ox * d_adj * f1;
    double y1 = endY - oy * d_adj * f1;
    double z1 = endZ - oz * d_adj * f1;

    double x2 = endX - ox * d_adj * f2;
    double y2 = endY - oy * d_adj * f2;
    double z2 = endZ - oz * d_adj * f2;

    QSMcylinder c;
    c.startX = x1;
    c.startY = y1;
    c.startZ = z1;
    c.endX   = x2;
    c.endY   = y2;
    c.endZ   = z2;
    c.cyl_ID     = next_id;
    c.parent_ID  = (i == 1 ? root->cyl_ID : prev_id);
    c.axis_ID    = 1;
    c.branch_order = 1;
    c.subtree_length = root_subtree + d_adj - actual_L * (nseg - i + 1);

    add_cylinder(c);

    prev_id = next_id;
    next_id--;
  }
}

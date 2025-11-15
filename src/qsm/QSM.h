#ifndef QSM_H
#define QSM_H

#include <unordered_map>
#include <vector>
#include <limits>
#include <cmath>

static constexpr double SUBTREE_LENGTH_UNSET    = -1.0;
static constexpr double SUBTREE_MAXZ_UNSET      = -1e300;
static constexpr double Z_EPS = 1e-9;

struct QSMcylinder
{
  double startX = 0;
  double startY = 0;
  double startZ = 0;
  double endX = 0;
  double endY = 0;
  double endZ = 0;
  double radius = 0;
  double subtree_length = SUBTREE_LENGTH_UNSET;
  double subtree_max_endZ = SUBTREE_MAXZ_UNSET;
  int cyl_ID = 0;
  int parent_ID = 0;
  int axis_ID = 0;
  int branch_order = 0;

  double length() const noexcept
  {
    double dx = endX - startX;
    double dy = endY - startY;
    double dz = endZ - startZ;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  double volume() const noexcept
  {
    double r2 = radius * radius;
    return M_PI * r2 * length();
  }
};

class QSM
{
public:
  QSM() = default;

  void build_from_cylinders(const std::vector<QSMcylinder>& cylinders);

  void compute_topology();
  void compute_architecture(int root_id = 1);

  const auto& cylinders() const { return cylinders_; }
  const auto& children_map() const { return children_map_; }

  size_t size() const { return cylinders_.size(); }

  // Range-based for loop supports
  auto begin() { return cylinders_.begin(); }
  auto end() { return cylinders_.end(); }
  auto begin() const { return cylinders_.begin(); }
  auto end() const { return cylinders_.end(); }

private:
  std::unordered_map<int, QSMcylinder> cylinders_;                 // cyl_ID -> cylinder
  std::unordered_map<int, std::vector<int>> children_map_;         // parent -> children

  // recursive helpers
  double compute_subtree_length(int node_id);
  double compute_subtree_max_z(int node_id);
  void assign_subtree_ids(int node_id, int current_axis_id, int current_branch_order, int &next_axis_id);
};


#endif // QSM_H

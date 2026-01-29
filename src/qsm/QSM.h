#ifndef QSM_H
#define QSM_H

#include <unordered_map>
#include <vector>
#include <limits>
#include <cmath>
#include <string>

#include "nanoflann/Adaptor.h"

using PointCloud = DataFrameAdaptor;

static constexpr double SUBTREE_LENGTH_UNSET    = -1.0;
static constexpr double SUBTREE_MAXZ_UNSET      = -1e300;
static constexpr double SUBTREE_VOLUME_UNSET    = -1;
static constexpr double RADIUS_UNSET            = -1.0;
static constexpr double Z_EPS = 1e-9;

struct QSMcylinder
{
  double startX = 0;
  double startY = 0;
  double startZ = 0;
  double endX = 0;
  double endY = 0;
  double endZ = 0;
  double radius = RADIUS_UNSET;
  double conic_allometry = RADIUS_UNSET;
  double subtree_length = SUBTREE_LENGTH_UNSET;
  double subtree_max_endZ = SUBTREE_MAXZ_UNSET;
  double subtree_volume = SUBTREE_VOLUME_UNSET;
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

  void add_cylinder(const QSMcylinder& c);

  void compute_topology();
  void smooth_skeleton(int niter, double th);
  void compute_architecture(int root_id = 1, bool use_volume = true);
  void prolongate(double d, double L = 0.1);
  void measure_radii(const PointCloud& tree, float sarc = 180, float sins = 0.2, float sinl = 0.3, float srmeas = 0.05);
  void polynomial_fitting(double tip_radius);
  void reconstruct_missing_radii(double tip_radius);
  void write(const std::string& filename, bool binary = true) const;

  const auto& cylinders() const { return cylinders_; }
  const auto& children_map() const { return children_map_; }

  size_t size() const { return cylinders_.size(); }

  // Range-based for loop supports
  auto begin() { return cylinders_.begin(); }
  auto end() { return cylinders_.end(); }
  auto begin() const { return cylinders_.begin(); }
  auto end() const { return cylinders_.end(); }

  QSMcylinder& get_cylinder_by_id(int cyl_id);

private:
  std::unordered_map<int, QSMcylinder> cylinders_;                 // cyl_ID -> cylinder
  std::unordered_map<int, std::vector<int>> children_map_;         // parent -> children

  std::vector<const QSMcylinder*> main_axis() const;

  // recursive helpers
  double compute_subtree_length(int node_id);
  double compute_subtree_max_z(int node_id);
  double compute_subtree_volume(int node_id);
  void assign_subtree_ids(int node_id, int current_axis_id, int current_branch_order, int &next_axis_id, bool use_volume);

  // mesh and write
  void build_mesh(std::vector<std::array<double,3>>& vertices,  std::vector<std::array<int,3>>& faces, int resolution = 16) const;
  void write_ply(const std::string& filename, bool binary) const;
  void write_stl(const std::string& filename, bool binary) const;
  void write_obj(const std::string& filename) const;
  void write_csv(const std::string& filename) const;

  // misc
  double conic_allometry(double tip_radius, double wi, double w0, double r0) const;
};


#endif // QSM_H

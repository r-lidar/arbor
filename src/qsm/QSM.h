#ifndef QSM_H
#define QSM_H

#include <unordered_map>
#include <vector>
#include <limits>
#include <cmath>
#include <string>

#include "Adaptor.h"

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

  inline double length() const noexcept
  {
    double dx = endX - startX;
    double dy = endY - startY;
    double dz = endZ - startZ;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  inline double volume() const noexcept
  {
    double r2 = radius * radius;
    return M_PI * r2 * length();
  }

  // Compute distance on the fly
  inline double distance(double px, double py, double pz) const
  {
    // 1. Compute Vector AB (Cylinder Axis)
    double ab_x = endX - startX;
    double ab_y = endY - startY;
    double ab_z = endZ - startZ;

    // 2. Compute Squared Length of AB
    double ab_sq = ab_x*ab_x + ab_y*ab_y + ab_z*ab_z;
    if(ab_sq < 1e-12) ab_sq = 1e-12; // Avoid division by zero

    // 3. Compute Vector AP (Start to Point)
    double ap_x = px - startX;
    double ap_y = py - startY;
    double ap_z = pz - startZ;

    // 4. Project AP onto AB (Dot Product)
    double t = (ap_x*ab_x + ap_y*ab_y + ap_z*ab_z) / ab_sq;

    // 5. Clamp t to [0, 1] to stay within the segment
    if (t < 0.0) t = 0.0;
    else if (t > 1.0) t = 1.0;

    // 6. Find Closest Point on Axis
    double c_x = startX + t * ab_x;
    double c_y = startY + t * ab_y;
    double c_z = startZ + t * ab_z;

    // 7. Distance from Point to Axis
    double dx = px - c_x;
    double dy = py - c_y;
    double dz = pz - c_z;

    double dist_axis = std::sqrt(dx*dx + dy*dy + dz*dz);

    // 8. Subtract radius for surface distance (clamped to 0)
    return (dist_axis > radius) ? (dist_axis - radius) : 0.0;
  }
};

struct Axe
{
  using container_type = std::vector<QSMcylinder*>;
  using iterator       = container_type::iterator;
  using const_iterator = container_type::const_iterator;

  container_type cylinders_;

  void add_cylinder(QSMcylinder* cyl) { cylinders_.push_back(cyl); }

  void sort()
  {
    // Sort cylinders from root to tip
    std::sort(cylinders_.begin(), cylinders_.end(),
      [](const QSMcylinder* a, const QSMcylinder* b) {  return a->subtree_length > b->subtree_length; });
  }

  bool need_reconstruction() const
  {
    for (const QSMcylinder* cyl : cylinders_)
    {
      if (cyl->radius == RADIUS_UNSET)
        return true;
    }
    return false;
  }

  void scale(double factor)
  {
    for (QSMcylinder* cyl : cylinders_) { cyl->radius *= factor; }
  }

  // --- vector-like access ---
  std::size_t size() const noexcept { return cylinders_.size(); }
  bool empty() const noexcept { return cylinders_.empty(); }
  QSMcylinder*& operator[](std::size_t i) { return cylinders_[i]; }
  const QSMcylinder* operator[](std::size_t i) const { return cylinders_[i]; }
  iterator begin() noexcept { return cylinders_.begin(); }
  iterator end()   noexcept { return cylinders_.end(); }
  const_iterator begin() const noexcept { return cylinders_.begin(); }
  const_iterator end()   const noexcept { return cylinders_.end(); }
  const_iterator cbegin() const noexcept { return cylinders_.cbegin(); }
  const_iterator cend()   const noexcept { return cylinders_.cend(); }
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
  void tmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,3>>& faces, int resolution = 16) const;
  void qmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,4>>& faces, int resolution = 16) const;
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
  void mesh(std::vector<std::array<double,3>>& vertices, int resolution) const;
  void write_ply(const std::string& filename, bool binary) const;
  void write_stl(const std::string& filename, bool binary) const;
  void write_obj(const std::string& filename) const;
  void write_csv(const std::string& filename) const;

  // misc
  double conic_allometry(double tip_radius, double wi, double w0, double r0) const;
};


#endif // QSM_H

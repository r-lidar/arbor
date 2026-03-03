#ifndef QSM_H
#define QSM_H

#include <unordered_map>
#include <vector>
#include <limits>
#include <cmath>
#include <string>

#include "PointCloud.h"
#include "QSMcylinder.h"

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

  static std::vector<std::pair<int, double>> layers(const PointCloud& points, double D);
  static std::vector<std::pair<int, double>> clusters(const PointCloud& points, const std::vector<std::pair<int, double>>&, double cl_dist);
  static PointCloud clean_tree_butt(const PointCloud&);

  void build_skeleton(const PointCloud&, const std::vector<std::pair<int, int>>& iter_cluster, double max_d);
  void compute_topology();
  void compute_architecture(int root_id = 1, bool use_volume = true);
  void smooth_skeleton(int niter, double th);
  void detect_weird_butt(double thresh = 50.0, int window = 4);
  void estimate_prolongation(const PointCloud& tree);
  void prolongate(double d, double L = 0.1);
  void construct_radii(const PointCloud& tree, double tip_radius = 0.0025);
  void measure_radii(const PointCloud& tree, float sarc = 180, float sins = 0.2, float sinl = 0.3, float srmeas = 0.05);
  void polynomial_fitting(double tip_radius = 0.0025);
  void reconstruct_missing_radii(double tip_radius);
  void tmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,3>>& faces, int resolution = 16) const;
  void qmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,4>>& faces, int resolution = 16) const;
  void write(const std::string& filename, bool binary = true) const;
  void conic_allometry(double R0, double tip_radius = 0.0025);
  void fix_multiple_root();
  void shift(double tx, double ty, double tz)
  {
    for (auto& kv : cylinders_)
    {
      kv.second.startX += tx;
      kv.second.endX += tx;
      kv.second.startY += ty;
      kv.second.endY += ty;
      kv.second.startZ += tz;
      kv.second.endZ += tz;
    }
  }

  const auto& cylinders() const { return cylinders_; }
  const auto& children_map() const { return children_map_; }

  size_t size() const { return cylinders_.size(); }

  // Range-based for loop supports
  auto begin() { return cylinders_.begin(); }
  auto end() { return cylinders_.end(); }
  auto begin() const { return cylinders_.begin(); }
  auto end() const { return cylinders_.end(); }

  QSMcylinder& get_cylinder_by_id(int cyl_id);

  double prolongation_distance = 0;

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
  void remove_disconnected_branches();
  double conic_allometry(double tip_radius, double wi, double w0, double r0) const;
};


#endif // QSM_H

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

class QSMbuilder;

class QSM
{
public:
  QSM() = default;

  // Accessor
  void add_cylinder(const QSMcylinder& c);
  const auto& cylinders() const { return cylinders_; }
  const auto& children_map() const { return children_map_; }
  std::vector<const QSMcylinder*> main_axis() const;
  size_t size() const { return cylinders_.size(); }

  // Mesh
  void tmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,3>>& faces, int resolution = 16) const;
  void qmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,4>>& faces, int resolution = 16) const;

  // I/O
  void write(const std::string& filename, bool binary = true) const;

  // Range-based for loop supports
  auto begin() { return cylinders_.begin(); }
  auto end() { return cylinders_.end(); }
  auto begin() const { return cylinders_.begin(); }
  auto end() const { return cylinders_.end(); }

  QSMcylinder& get_cylinder_by_id(int cyl_id);

private:
  friend class QSMbuilder;
  std::unordered_map<int, QSMcylinder> cylinders_;                 // cyl_ID -> cylinder
  std::unordered_map<int, std::vector<int>> children_map_;         // parent -> children

  // Mesh
  void mesh(std::vector<std::array<double,3>>& vertices, int resolution) const;

  // I/O
  void write_ply(const std::string& filename, bool binary) const;
  void write_stl(const std::string& filename, bool binary) const;
  void write_obj(const std::string& filename) const;
  void write_csv(const std::string& filename) const;
};


#endif // QSM_H

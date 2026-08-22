/**
 * @file PointCloudDefault.h
 * Project: Arbor
 *
 * Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#ifndef POINTCLOUDDEFAULT_H
#define POINTCLOUDDEFAULT_H

#ifndef USING_R

#include "PointCloudBase.h"
#include <array>
#include <span>

#define POINT_CLOUD_ATTR(interface_type, internal_type, name, container, has_func, error_msg)  \
inline interface_type get_##name(const size_t idx) const {                                     \
  if (!has_func()) throw std::runtime_error(error_msg);                                        \
  return static_cast<interface_type>(container[idx]);                                          \
}                                                                                              \
inline void set_##name(const size_t idx, interface_type v) {                                   \
  if (!has_func()) throw std::runtime_error(error_msg);                                        \
  container[idx] = static_cast<internal_type>(v);                                              \
}

class PointCloudDefault : public BasePointCloud<PointCloudDefault>
{
public:
  using BasePointCloud<PointCloudDefault>::subset;

  // Constructors / destructor
  PointCloudDefault();
  PointCloudDefault(size_t n, bool init_attributes = false);
  PointCloudDefault(const PointCloudDefault&) = default;
  PointCloudDefault(PointCloudDefault&&) = default;
  PointCloudDefault& operator=(const PointCloudDefault&) = default;
  PointCloudDefault& operator=(PointCloudDefault&&) = default;
  ~PointCloudDefault() = default;

  // Insertion
  void add_point(float x, float y, float z)
  {
    n_points++;
    true_n_points++;
    coords.push_back({x,y,z});
    treeid.push_back(0);
    foliage.push_back(1);
    classif.push_back(0);
    passage.push_back(0);
    hag.push_back(0);
    pwood.push_back(0);
    rgb.push_back({200, 200, 200});
    kdtree.reset();
  }

  // Merge
  PointCloudDefault& operator+=(const PointCloudDefault& other);
  PointCloudDefault  operator+(const PointCloudDefault& other) const;

  // --- Nanoflann KD-tree interface ---
  inline size_t kdtree_get_point_count() const { return n_points; }
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const
  {
    switch (dim)
    {
    case 0: return static_cast<double>(coords[idx].x);
    case 1: return static_cast<double>(coords[idx].y);
    case 2: return static_cast<double>(coords[idx].z);
    default: throw std::runtime_error("Invalid dimension");
    }
  }
  template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }

  // --- Num. points ---
  inline size_t point_count() const { return n_points; }
  inline size_t size()        const { return n_points; }
  inline size_t true_size()   const { return true_n_points; }

  // --- Geometry access ---
  inline void get_point(const size_t idx, double* q) const
  {
    q[0] = static_cast<double>(coords[idx].x);
    q[1] = static_cast<double>(coords[idx].y);
    q[2] = static_cast<double>(coords[idx].z);
  }

  inline double get_x(const size_t idx) const { return static_cast<double>(coords[idx].x); }
  inline double get_y(const size_t idx) const { return static_cast<double>(coords[idx].y); }
  inline double get_z(const size_t idx) const { return static_cast<double>(coords[idx].z); }
  inline void   set_x(const size_t idx, double v) { coords[idx].x = static_cast<float>(v); }
  inline void   set_y(const size_t idx, double v) { coords[idx].y = static_cast<float>(v); }
  inline void   set_z(const size_t idx, double v) { coords[idx].z = static_cast<float>(v); }

  // --- Optional attribute access ---
  inline bool has_red()      const { return !rgb.empty(); }
  inline bool has_blue()     const { return !rgb.empty(); }
  inline bool has_green()    const { return !rgb.empty(); }
  inline bool has_rgb()      const { return !rgb.empty(); }
  inline bool has_hag()      const { return !hag.empty(); }
  inline bool has_treeid()   const { return !treeid.empty(); }
  inline bool has_pwood()    const { return !pwood.empty(); }
  inline bool has_foliage()  const { return !foliage.empty(); }
  inline bool has_passage()  const { return !passage.empty(); }
  inline bool has_class()    const { return !classif.empty(); }
  inline bool has_userdata() const { return !userdata.empty(); }

  POINT_CLOUD_ATTR(int,    int,      treeid,         treeid,   has_treeid,   "Instance segmentation not available")
  POINT_CLOUD_ATTR(double, double,   pwood,          pwood,    has_pwood,    "Wood likelihood data not available")
  POINT_CLOUD_ATTR(int,    uint8_t,  foliage,        foliage,  has_foliage,  "Semantic segmentation not available")
  POINT_CLOUD_ATTR(double, float,    hag,            hag,      has_hag,      "HAG data not available")
  POINT_CLOUD_ATTR(int,    int,      passage,        passage,  has_passage,  "Passage not available")
  POINT_CLOUD_ATTR(int,    uint16_t, classification, classif,  has_class,    "Classification data not available")
  POINT_CLOUD_ATTR(int,    uint8_t,  userdata,       userdata, has_userdata, "Classification data not available")

  #undef POINT_CLOUD_ATTR

  inline void set_color(size_t i, uint8_t r, uint8_t g, uint8_t b) { rgb[i] = {r, g, b}; }
  inline std::array<uint8_t, 3> get_color(size_t i) const
  {
    const RGB& color = rgb[i];
    return {color.r, color.g, color.b};
  }

  inline bool is_wood(const size_t idx) const {
    return get_foliage(idx) == 0;
  }

  inline bool is_ground(const size_t idx) const {
    return get_classification(idx) == 2;
  }

  // --- In-place transforms ---
  void translate(double x, double y, double z) ;
  void scale(double x, double y, double z) ;

  // --- Subset ---
  PointCloudDefault subset(const std::vector<unsigned int>& indices, bool xyz_only = false) const;

  // --- View ----
  std::span<float> coord_view()
  {
    if (coords.empty()) return {};
    return { reinterpret_cast<float*>(coords.data()), coords.size() * 3};
  }

  std::span<uint8_t> rgb_view()
  {
    if (rgb.empty()) return {};
    return { reinterpret_cast<uint8_t*>(rgb.data()), rgb.size() * 3};
  }

private:
  void cleanup();
  void init();
  void safe_alloc(size_t n, bool alloc_attrs);
  void swap_points(size_t i, size_t j)
  {
    std::swap(coords[i], coords[j]);
    if (has_rgb())     std::swap(rgb[i], rgb[j]);
    if (has_treeid())  std::swap(treeid[i], treeid[j]);
    if (has_foliage()) std::swap(foliage[i], foliage[j]);
    if (has_passage()) std::swap(passage[i], passage[j]);
    if (has_class())   std::swap(classif[i], classif[j]);
    if (has_hag())     std::swap(hag[i], hag[j]);
    if (has_pwood())   std::swap(pwood[i], pwood[j]);
  }

private:
  struct Vec3 { float x,y,z; };
  struct RGB { uint8_t r,g,b; };
  std::vector<Vec3> coords;
  std::vector<RGB> rgb;
  std::vector<int> treeid;
  std::vector<uint8_t> foliage;
  std::vector<uint16_t> classif;
  std::vector<int> passage;
  std::vector<float> hag;
  std::vector<float> pwood;
  std::vector<uint8_t> userdata;
};

#endif

#endif

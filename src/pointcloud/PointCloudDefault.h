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

#define POINT_CLOUD_ATTR(interface_type, name, container, flag, error_msg)                       \
inline interface_type get_##name(const size_t idx) const {                                     \
  if (!has(flag)) throw std::runtime_error(error_msg);                                         \
  return static_cast<interface_type>(container[idx]);                                          \
}                                                                                              \
inline void set_##name(const size_t idx, interface_type v) {                                   \
  if (!has(flag)) throw std::runtime_error(error_msg);                                         \
  container[idx] = static_cast<typename decltype(container)::value_type>(v);                   \
}

class PointCloudDefault : public BasePointCloud<PointCloudDefault>
{
public:
  using BasePointCloud<PointCloudDefault>::subset;

  // Constructors / destructor
  // XYZ is always allocated. attrs selects the optional attributes, filled with their default value.
  PointCloudDefault();
  explicit PointCloudDefault(size_t n, Attr attrs = Attr::None);
  PointCloudDefault(const PointCloudDefault&) = default;
  PointCloudDefault(PointCloudDefault&&) = default;
  PointCloudDefault& operator=(const PointCloudDefault&) = default;
  PointCloudDefault& operator=(PointCloudDefault&&) = default;
  ~PointCloudDefault() = default;

  // Insertion. The present attributes get their default value.
  void add_point(float x, float y, float z);

  // --- Attribute management ---
  inline Attr attributes() const { return present; }
  inline bool has(Attr a) const { return (present & a) == a; }

  // Allocates the attributes of a that are missing, filled with their default value. Inside a ScopedPartition the
  // hidden points are allocated too.
  void allocate(Attr a);
  // Frees the attributes of a.
  void release(Attr a);

  // Merge. operator+= throws if the two point clouds do not have the same attributes. merge() keeps only the
  // attributes of keep (the others are released) and throws if one of the point clouds lacks one of them.
  // Not inside a ScopedPartition.
  PointCloudDefault& merge(const PointCloudDefault& other, Attr keep);
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
  inline bool has_red()      const { return has(Attr::RGB); }
  inline bool has_blue()     const { return has(Attr::RGB); }
  inline bool has_green()    const { return has(Attr::RGB); }
  inline bool has_rgb()      const { return has(Attr::RGB); }
  inline bool has_intensity() const { return has(Attr::Intensity); }
  inline bool has_gpstime()  const { return has(Attr::GPStime); }
  inline bool has_return_number() const { return has(Attr::ReturnNumber); }
  inline bool has_number_of_returns() const { return has(Attr::NumberOfReturns); }
  inline bool has_hag()      const { return has(Attr::HAG); }
  inline bool has_treeid()   const { return has(Attr::TreeID); }
  inline bool has_pwood()    const { return has(Attr::PWood); }
  inline bool has_foliage()  const { return has(Attr::Foliage); }
  inline bool has_passage()  const { return has(Attr::Passage); }
  inline bool has_class()    const { return has(Attr::Classification); }
  inline bool has_userdata() const { return has(Attr::UserData); }

  POINT_CLOUD_ATTR(int,    treeid,            treeid,            Attr::TreeID,          "Instance segmentation not available")
  POINT_CLOUD_ATTR(double, pwood,             pwood,             Attr::PWood,           "Wood likelihood data not available")
  POINT_CLOUD_ATTR(int,    foliage,           foliage,           Attr::Foliage,         "Semantic segmentation not available")
  POINT_CLOUD_ATTR(double, hag,               hag,               Attr::HAG,             "HAG data not available")
  POINT_CLOUD_ATTR(int,    passage,           passage,           Attr::Passage,         "Passage not available")
  POINT_CLOUD_ATTR(int,    classification,    classif,           Attr::Classification,  "Classification data not available")
  POINT_CLOUD_ATTR(int,    userdata,          userdata,          Attr::UserData,        "UserData not available")
  POINT_CLOUD_ATTR(int,    intensity,         intensity,         Attr::Intensity,       "Intensity not available")
  POINT_CLOUD_ATTR(double, gpstime,           gpstime,           Attr::GPStime,         "gpstime not available")
  POINT_CLOUD_ATTR(int,    return_number,     return_number,     Attr::ReturnNumber,    "ReturnNumber not available")
  POINT_CLOUD_ATTR(int,    number_of_returns, number_of_returns, Attr::NumberOfReturns, "NumberOfReturns not available")

  #undef POINT_CLOUD_ATTR

  inline void set_color(size_t i, uint8_t r, uint8_t g, uint8_t b)
  {
    if (!has_rgb()) throw std::runtime_error("RGB not available");
    rgb[i] = {r, g, b};
  }
  inline std::array<uint8_t, 3> get_color(size_t i) const
  {
    if (!has_rgb()) throw std::runtime_error("RGB not available");
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
  // Copies the attributes of keep that are present.
  PointCloudDefault subset(const std::vector<unsigned int>& indices, Attr keep = Attr::All) const;

  // --- Filter ---
  // Unlike subset(), which returns a new point cloud, removes in place the points whose `keep` is false (one value
  // per point). The kept points keep their order and every attribute; the memory of the removed points is
  // released, one attribute at a time. Not inside a ScopedPartition.
  void filter(const std::vector<bool>& keep);

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
  void swap_points(size_t i, size_t j) override;

  struct Vec3 { float x,y,z; };
  struct RGB { uint8_t r,g,b; };

  // The table of the optional attributes: calls f(flag, member pointer, default value) for each of them.
  // Every operation that touches all the attributes (allocation, subset, merge, filter, swap) loops over it.
  template <class F>
  static void visit_attributes(F&& f)
  {
    f(Attr::RGB,             &PointCloudDefault::rgb,               RGB{128, 128, 128});
    f(Attr::Intensity,       &PointCloudDefault::intensity,         uint16_t{0});
    f(Attr::GPStime,         &PointCloudDefault::gpstime,           0.0);
    f(Attr::ReturnNumber,    &PointCloudDefault::return_number,     uint8_t{1});
    f(Attr::NumberOfReturns, &PointCloudDefault::number_of_returns, uint8_t{1});
    f(Attr::Classification,  &PointCloudDefault::classif,           uint16_t{0});
    f(Attr::UserData,        &PointCloudDefault::userdata,          uint8_t{0});
    f(Attr::HAG,             &PointCloudDefault::hag,               0.0f);
    f(Attr::PWood,           &PointCloudDefault::pwood,             0.0f);
    f(Attr::Foliage,         &PointCloudDefault::foliage,           uint8_t{1});
    f(Attr::Passage,         &PointCloudDefault::passage,           0);
    f(Attr::TreeID,          &PointCloudDefault::treeid,            0);
  }

private:
  Attr present = Attr::None;
  std::vector<Vec3> coords;
  std::vector<RGB> rgb;
  std::vector<uint16_t> intensity;
  std::vector<double> gpstime;
  std::vector<uint8_t> return_number;
  std::vector<uint8_t> number_of_returns;
  std::vector<uint16_t> classif;
  std::vector<uint8_t> userdata;
  std::vector<float> hag;
  std::vector<float> pwood;
  std::vector<uint8_t> foliage;
  std::vector<int> passage;
  std::vector<int> treeid;
};

#endif

#endif

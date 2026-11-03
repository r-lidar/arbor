/**
 * @file PointCloudDataframe.h
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

#ifndef POINTCLOUDDATAFRAME_H
#define POINTCLOUDDATAFRAME_H

#include "PointCloudBase.h"
#include <Rcpp.h>

#ifdef USING_R

#define POINT_CLOUD_ATTR(type, name, container, has_func, error_msg) \
inline type get_##name(const size_t idx) const {                     \
  if (!has_func()) throw std::runtime_error(error_msg);              \
  return container[idx];                                             \
}                                                                    \
inline void set_##name(const size_t idx, type v) {                   \
  if (!has_func()) throw std::runtime_error(error_msg);              \
  container[idx] = v;                                                \
}                                                                    \

class PointCloudDataFrame : public BasePointCloud<PointCloudDataFrame>
{
public:
  using BasePointCloud<PointCloudDataFrame>::subset;

  // Constructors / destructor
  // XYZ is always allocated. attrs selects the optional attributes, filled with their default value. The set of
  // attributes is fixed at construction: there is no allocate().
  PointCloudDataFrame();
  explicit PointCloudDataFrame(size_t n, Attr attrs = Attr::None);

  // Rule of five
  PointCloudDataFrame(const PointCloudDataFrame& other);
  PointCloudDataFrame(PointCloudDataFrame&& other) noexcept;
  PointCloudDataFrame& operator=(const PointCloudDataFrame other) noexcept;
  ~PointCloudDataFrame() ;

  // --- Attribute management ---
  inline Attr attributes() const { return present; }
  inline bool has(Attr a) const { return (present & a) == a; }

  // Not supported: the attributes are either the columns of the R data.frame or allocated by the constructor.
  [[noreturn]] void allocate(Attr a);
  [[noreturn]] void release(Attr a);

  // Merge (owned memory only). operator+= throws if the two point clouds do not have the same attributes. merge()
  // keeps only the attributes of keep (the others are freed) and throws if one of the point clouds lacks one of them.
  PointCloudDataFrame& merge(const PointCloudDataFrame& other, Attr keep);
  PointCloudDataFrame& operator+=(const PointCloudDataFrame& other);
  PointCloudDataFrame  operator+(const PointCloudDataFrame& other) const;

  // --- Nanoflann KD-tree interface ---
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const { return coords[dim][idx]; }

  // --- Geometry access ---
  inline void get_point(const size_t idx, double* q) const { for (size_t d = 0; d < 3; ++d) q[d] = coords[d][idx];}

  inline double get_x(const size_t idx) const { return coords[0][idx]; }
  inline double get_y(const size_t idx) const { return coords[1][idx]; }
  inline double get_z(const size_t idx) const { return coords[2][idx]; }
  inline void   set_x(const size_t idx, double v) { coords[0][idx] = v; }
  inline void   set_y(const size_t idx, double v) { coords[1][idx] = v; }
  inline void   set_z(const size_t idx, double v) { coords[2][idx] = v; }

  // --- Optional attribute access ---
  inline bool has_hag()      const { return has(Attr::HAG); }
  inline bool has_treeid()   const { return has(Attr::TreeID); }
  inline bool has_pwood()    const { return has(Attr::PWood); }
  inline bool has_foliage()  const { return has(Attr::Foliage); }
  inline bool has_passage()  const { return has(Attr::Passage); }
  inline bool has_class()    const { return has(Attr::Classification); }
  inline bool has_red()      const { return has(Attr::RGB); }
  inline bool has_green()    const { return has(Attr::RGB); }
  inline bool has_blue()     const { return has(Attr::RGB); }
  inline bool has_rgb()      const { return has(Attr::RGB); }
  inline bool has_userdata() const { return has(Attr::UserData); }
  inline bool has_intensity() const { return has(Attr::Intensity); }
  inline bool has_gpstime()  const { return has(Attr::GPStime); }
  inline bool has_return_number() const { return has(Attr::ReturnNumber); }
  inline bool has_number_of_returns() const { return has(Attr::NumberOfReturns); }

  POINT_CLOUD_ATTR(int, treeid, treeid, has_treeid, "Instance segmentation not available")
  POINT_CLOUD_ATTR(double, pwood, pwood, has_pwood, "Wood likelihood data not available")
  POINT_CLOUD_ATTR(int, foliage, foliage, has_foliage, "Semantic segmentation not available")
  POINT_CLOUD_ATTR(double, hag, hag, has_hag, "HAG not available")
  POINT_CLOUD_ATTR(int, passage, passage, has_passage, "Passage not available")
  POINT_CLOUD_ATTR(int, classification, classif, has_class, "Classification not available")
  POINT_CLOUD_ATTR(int, red, red, has_red, "RGB not available")
  POINT_CLOUD_ATTR(int, green, green, has_green, "RGB not available")
  POINT_CLOUD_ATTR(int, blue, blue, has_blue, "RGB not available")
  POINT_CLOUD_ATTR(int, userdata, userdata, has_userdata, "UserData not available")
  POINT_CLOUD_ATTR(int, intensity, intensity, has_intensity, "Intensity not available")
  POINT_CLOUD_ATTR(double, gpstime, gpstime, has_gpstime, "gpstime not available")
  POINT_CLOUD_ATTR(int, return_number, return_number, has_return_number, "ReturnNumber not available")
  POINT_CLOUD_ATTR(int, number_of_returns, number_of_returns, has_number_of_returns, "NumberOfReturns not available")

  #undef POINT_CLOUD_ATTR

  inline void   set_color(size_t i, uint8_t r, uint8_t g, uint8_t b)
  {
    set_red(i,   static_cast<int>(r) * 255);
    set_green(i, static_cast<int>(g) * 255);
    set_blue(i,  static_cast<int>(b) * 255);
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
  PointCloudDataFrame subset(const std::vector<unsigned int>&  indices, Attr keep = Attr::All) const;

  // --- Used only in R ----
  explicit PointCloudDataFrame(const Rcpp::DataFrame& df);
  bool linked_to_dataframe() const { return !owns_memory; }

private:
  void cleanup();
  void swap(PointCloudDataFrame& first, PointCloudDataFrame& second) noexcept;
  void safe_alloc(size_t n, Attr attrs);
  void swap_points(size_t i, size_t j) override;

  // The table of the optional columns: calls f(flag, member pointer, default value, R column name) for each of
  // them. RGB is three columns sharing the same flag. Every operation that touches all the attributes (allocation,
  // copy, subset, merge, swap, cleanup) loops over it. The type of the member gives the R type: int* is INTSXP,
  // double* is REALSXP.
  template <class F>
  static void visit_attributes(F&& f)
  {
    f(Attr::RGB,             &PointCloudDataFrame::red,               128 * 255, "R");
    f(Attr::RGB,             &PointCloudDataFrame::green,             128 * 255, "G");
    f(Attr::RGB,             &PointCloudDataFrame::blue,              128 * 255, "B");
    f(Attr::Intensity,       &PointCloudDataFrame::intensity,         0,   "Intensity");
    f(Attr::GPStime,         &PointCloudDataFrame::gpstime,           0.0, "gpstime");
    f(Attr::ReturnNumber,    &PointCloudDataFrame::return_number,     1,   "ReturnNumber");
    f(Attr::NumberOfReturns, &PointCloudDataFrame::number_of_returns, 1,   "NumberOfReturns");
    f(Attr::Classification,  &PointCloudDataFrame::classif,           0,   "Classification");
    f(Attr::UserData,        &PointCloudDataFrame::userdata,          0,   "UserData");
    f(Attr::HAG,             &PointCloudDataFrame::hag,               0.0, "hag");
    f(Attr::PWood,           &PointCloudDataFrame::pwood,             0.0, "pwood");
    f(Attr::Foliage,         &PointCloudDataFrame::foliage,           1,   "foliage");
    f(Attr::Passage,         &PointCloudDataFrame::passage,           0,   "passage");
    f(Attr::TreeID,          &PointCloudDataFrame::treeid,            0,   "treeID");
  }

private:
  double* coords[3] = {nullptr, nullptr, nullptr};

  // Optional attributes (nullptr if absent)
  int*    red               = nullptr;
  int*    green             = nullptr;
  int*    blue              = nullptr;
  int*    intensity         = nullptr;
  double* gpstime           = nullptr;
  int*    return_number     = nullptr;
  int*    number_of_returns = nullptr;
  int*    classif           = nullptr;
  int*    userdata          = nullptr;
  double* hag               = nullptr;
  double* pwood             = nullptr;
  int*    foliage           = nullptr;
  int*    passage           = nullptr;
  int*    treeid            = nullptr;

  Attr present = Attr::None;
  bool owns_memory = true;
  struct RawColumn { void* ptr; SEXPTYPE type; };
  std::vector<RawColumn> generic_columns;
};

#endif

#endif

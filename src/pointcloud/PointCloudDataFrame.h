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
  PointCloudDataFrame();
  PointCloudDataFrame(size_t n, bool init_attributes = false);

  // Rule of five
  PointCloudDataFrame(const PointCloudDataFrame& other);
  PointCloudDataFrame(PointCloudDataFrame&& other) noexcept;
  PointCloudDataFrame& operator=(const PointCloudDataFrame other) noexcept;
  ~PointCloudDataFrame() ;

  // Merge
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
  inline bool has_hag()     const { return hag != nullptr; }
  inline bool has_treeid()  const { return treeid != nullptr; }
  inline bool has_pwood()   const { return pwood != nullptr; }
  inline bool has_foliage() const { return foliage != nullptr; }
  inline bool has_passage() const { return passage != nullptr; }
  inline bool has_class()   const { return classif != nullptr; }
  inline bool has_red()     const { return red != nullptr; }
  inline bool has_green()   const { return green != nullptr; }
  inline bool has_blue()    const { return blue != nullptr; }

  POINT_CLOUD_ATTR(int, treeid, treeid, has_treeid, "Instance segmentation not available")
  POINT_CLOUD_ATTR(double, pwood, pwood, has_pwood, "Wood likelihood data not available")
  POINT_CLOUD_ATTR(int, foliage, foliage, has_foliage, "Semantic segmentation not available")
  POINT_CLOUD_ATTR(double, hag, hag, has_hag, "HAG not available")
  POINT_CLOUD_ATTR(int, passage, passage, has_passage, "Passage not available")
  POINT_CLOUD_ATTR(int, classification, classif, has_class, "Classification not available")
  POINT_CLOUD_ATTR(int, red, red, has_red, "RGB not available")
  POINT_CLOUD_ATTR(int, green, green, has_green, "RGB not available")
  POINT_CLOUD_ATTR(int, blue, blue, has_blue, "RGB not available")

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
  PointCloudDataFrame subset(const std::vector<unsigned int>&  indices, bool xyz_only = false) const;

  // --- Used only in R ----
  explicit PointCloudDataFrame(const Rcpp::DataFrame& df);
  bool linked_to_dataframe() const { return !owns_memory; }

private:
  void cleanup();
  void swap(PointCloudDataFrame& first, PointCloudDataFrame& second) noexcept;
  void init();
  void safe_alloc(size_t n, bool alloc_attrs);
  void swap_points(size_t i, size_t j) override
  {
    for (int d = 0; d < 3; ++d) std::swap(coords[d][i], coords[d][j]);
    if (treeid)  std::swap(treeid[i], treeid[j]);
    if (foliage) std::swap(foliage[i], foliage[j]);
    if (passage) std::swap(passage[i], passage[j]);
    if (classif) std::swap(classif[i], classif[j]);
    if (hag)     std::swap(hag[i], hag[j]);
    if (pwood)   std::swap(pwood[i], pwood[j]);
    if (red)     std::swap(red[i],   red[j]);
    if (green)   std::swap(green[i], green[j]);
    if (blue)    std::swap(blue[i],  blue[j]);
  }
  template <typename T>
  void merge_attribute(T*& current, const T* other, size_t old_n, size_t other_n, bool other_has_it)
  {
    if (current && other_has_it)
    {
      // Both have the attribute: Reallocate and merge
      T* next = new T[old_n + other_n];
      std::memcpy(next, current, old_n * sizeof(T));
      std::memcpy(next + old_n, other, other_n * sizeof(T));

      delete[] current;
      current = next;
    }
    else if (current)
    {
      // Only 'this' has it: Drop the attribute to maintain consistency
      delete[] current;
      current = nullptr;
    }
  }

private:
  double* coords[3] = {nullptr, nullptr, nullptr};

  // Optional attributes (nullptr if absent)
  int*    treeid  = nullptr;
  int*    foliage = nullptr;
  int*    passage = nullptr;
  int*    classif = nullptr;
  int*    red     = nullptr;
  int*    green   = nullptr;
  int*    blue    = nullptr;
  double* hag     = nullptr;
  double* pwood   = nullptr;

  bool owns_memory = false;
};

#endif

#endif

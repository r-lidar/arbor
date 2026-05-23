/**
 * @file PointCloud.h
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

#ifndef ADAPTATOR_H
#define ADAPTATOR_H

#include <memory>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdint>

#include "nanoflann.h"

#ifdef USING_R
#include <Rcpp.h>
#else
#include <span>
#endif

class BasePointCloud
{
public:
  // For nanoflann
  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const { return false; }
  inline size_t kdtree_get_point_count() const { return size(); }

  // Generic
  inline size_t size()        const { return n_points; }
  inline size_t true_size()   const { return true_n_points; }

  // Partition like std::partition
  template <typename Predicate>
  size_t partition(Predicate pred)
  {
    if (n_points == 0) return 0;

    int64_t left = 0;
    int64_t right = static_cast<int64_t>(n_points) - 1;

    while (left <= right)
    {
      // Move left pointer until we find a point that should be at the back
      while (left <= right && pred(left))
      {
        left++;
      }
      // Move right pointer until we find a point that should be at the front
      while (left <= right && !pred(right))
      {
        right--;
      }

      if (left <= right)
      {
        swap_points(left, right);
        left++;
        right--;
      }
    }

    // Update n_points so nanoflann only sees the front section
    // Or keep n_points same and return 'left' as the new logical size
    size_t new_size = static_cast<size_t>(left);
    n_points = new_size;
    return new_size;
  }

protected:
  virtual void swap_points(size_t i, size_t j) = 0; // Needed for partition
  size_t n_points;
  size_t true_n_points;
};

/**
 * @brief Blueprint for a custom PointCloud implementation.
 * Users can copy-paste this and fill in their specific memory logic.
 */
/*class MinimalPointCloud : public BasePointCloud
{
public:
  // Constructors / destructor
  MinimalPointCloud() = default;
  MinimalPointCloud(size_t n, bool init_attributes = false);

  // Rule of five
  MinimalPointCloud(const MinimalPointCloud& other);               // Copy constructor
  MinimalPointCloud(MinimalPointCloud&& other) noexcept;           // Move constructor
  MinimalPointCloud& operator=(MinimalPointCloud other) noexcept;  // Assignment Operator (using Copy-and-Swap idiom)
  ~MinimalPointCloud() ;                                           // destructor

  // Merging operator
  MinimalPointCloud& operator+=(const MinimalPointCloud& other);
  MinimalPointCloud  operator+(const MinimalPointCloud& other) const;


  // Nanoflann KD-tree interface
  inline size_t kdtree_get_point_count() const { return 0; }
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const { return 0; }
  template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }

  // Num. points
  inline size_t size()        const { return 0; }
  inline size_t true_size()   const { return 0; }

  // Geometry access
  inline void   get_point(const size_t idx, double* q) const { }
  inline double get_x(const size_t idx) const { return 0; }
  inline double get_y(const size_t idx) const { return 0; }
  inline double get_z(const size_t idx) const { return 0; }
  inline void   set_x(const size_t idx, double v) { }
  inline void   set_y(const size_t idx, double v) { }
  inline void   set_z(const size_t idx, double v) { }

  // Attribute access. Each point must have an attribute hag/treeid/pwood/foliage/passage
  inline bool has_hag()     const { return true; }
  inline bool has_treeid()  const { return true; }
  inline bool has_pwood()   const { return true; }
  inline bool has_foliage() const { return true; }
  inline bool has_passage() const { return true; }

  inline int    get_treeid (const size_t idx) const { return 0; }
  inline void   set_treeid (const size_t idx, int v) { }
  inline double get_pwood  (const size_t idx) const { return 0; }
  inline void   set_pwood  (const size_t idx, double v) { }
  inline int    get_foliage(const size_t idx) const { return 0; }
  inline void   set_foliage(const size_t idx, int v) { }
  inline double get_hag    (const size_t idx) const { return 0; }
  inline void   set_hag    (const size_t idx, double v) { }
  inline int    get_passage(const size_t idx) const { return 0; }
  inline void   set_passage(const size_t idx, int v) { }
  inline void   set_ground(const size_t idx, bool v) { };
  inline bool   is_ground(const size_t idx) { return false; }

  // 0 wood 1 2 foliage
  inline bool is_wood(const size_t idx) const { return get_foliage(idx) == 0; }

  // In-place transforms. No need to be inlined
  void translate(double x, double y, double z);
  void scale(double x, double y, double z);

  // Subset
  MinimalPointCloud subset(const std::vector<bool>& keep, bool xyz_only = false) const;


private:
  // Custom memory layout
  void swap_points(size_t i, size_t j) { return; }
};*/

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

class PointCloudDataFrame : public BasePointCloud
{
  using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloudDataFrame>, PointCloudDataFrame, 3>;

public:
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
  PointCloudDataFrame subset(const std::vector<int>&  indices, bool xyz_only = false) const;
  PointCloudDataFrame subset(const std::vector<bool>& mask,    bool xyz_only = false) const;

  // --- Color ---
  void colorize_trees(bool darken_foliage = false);

  // --- Used only in R ----
  explicit PointCloudDataFrame(const Rcpp::DataFrame& df);
  bool linked_to_dataframe() const { return !owns_memory; }

  // --- Spatial queries -------
  void build_index() const;
  void knn(const double* query, int k, std::vector<unsigned int>&idx, std::vector<double>& sqdist) const;
  void radius_search(const double* query, double radius, std::vector<unsigned int>& indices, std::vector<double>& sqdist) const;


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

  mutable std::unique_ptr<KDTree> kdtree;
};

#else

#define POINT_CLOUD_ATTR(interface_type, internal_type, name, container, has_func, error_msg)  \
inline interface_type get_##name(const size_t idx) const {                                     \
  if (!has_func()) throw std::runtime_error(error_msg);                                        \
  return static_cast<interface_type>(container[idx]);                                          \
}                                                                                              \
inline void set_##name(const size_t idx, interface_type v) {                                   \
  if (!has_func()) throw std::runtime_error(error_msg);                                        \
  container[idx] = static_cast<internal_type>(v);                                              \
}

class PointCloudDefault : public BasePointCloud
{
public:
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
  inline bool has_rgb()     const { return !rgb.empty(); }
  inline bool has_hag()     const { return !hag.empty(); }
  inline bool has_treeid()  const { return !treeid.empty(); }
  inline bool has_pwood()   const { return !pwood.empty(); }
  inline bool has_foliage() const { return !foliage.empty(); }
  inline bool has_passage() const { return !passage.empty(); }
  inline bool has_class()   const { return !classif.empty(); }

  POINT_CLOUD_ATTR(int,    int,      treeid,         treeid,  has_treeid,  "Instance segmentation not available")
  POINT_CLOUD_ATTR(double, double,   pwood,          pwood,   has_pwood,   "Wood likelihood data not available")
  POINT_CLOUD_ATTR(int,    uint8_t,  foliage,        foliage, has_foliage, "Semantic segmentation not available")
  POINT_CLOUD_ATTR(double, float,    hag,            hag,     has_hag,     "HAG data not available")
  POINT_CLOUD_ATTR(int,    int,      passage,        passage, has_passage, "Passage not available")
  POINT_CLOUD_ATTR(int,    uint16_t, classification, classif, has_class,   "Classification data not available")

  #undef POINT_CLOUD_ATTR

  inline bool is_wood(const size_t idx) const {
      return get_foliage(idx) == 0;
  }

  inline bool is_ground(const size_t idx) const {
      return get_classification(idx) == 2;
  }

  // --- In-place transforms ---
  void translate(double x, double y, double z) ;
  void scale(double x, double y, double z) ;

  // --- Search ---
  std::vector<unsigned int> find_points_by_treeid(int tid) const;

  // --- Subset ---
  PointCloudDefault subset(const std::vector<int>& indices, bool xyz_only = false) const;
  PointCloudDefault subset(const std::vector<bool>& mask, bool xyz_only = false) const;
  PointCloudDefault subset_by_treeid(int tid) const;

  // --- Color ---
  void colorize_trees(bool darken_foliage = false);

  // --- Spatial queries -------
  void build_index() const;
  void knn(const double* query, int k, std::vector<unsigned int>&idx, std::vector<double>& sqdist) const;
  void radius_search(const double* query, double radius, std::vector<unsigned int>& indices, std::vector<double>& sqdist) const;

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
};
#endif

#ifdef USING_R
using PointCloud = PointCloudDataFrame;
#else
using PointCloud = PointCloudDefault;
#endif

#endif

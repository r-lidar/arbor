#ifndef ADAPTATOR_H
#define ADAPTATOR_H

#include <vector>
#include <string>
#include <stdexcept>

#ifdef USING_R
#include <Rcpp.h>
#endif

/**
 * @brief Blueprint for a custom PointCloud implementation.
 * Users can copy-paste this and fill in their specific memory logic.
 */
class MinimalPointCloud
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
  inline size_t point_count() const { return 0; }
  inline size_t size()        const { return 0; }

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

  // 0 wood 1 2 foliage
  inline bool is_wood(const size_t idx) const { return get_foliage(idx) == 0; }

  // In-place transforms. No need to be inlined
  void translate(double x, double y, double z);
  void scale(double x, double y, double z);

  // Subset
  MinimalPointCloud subset(const std::vector<bool>& keep, bool xyz_only = false) const;


private:
  // Custom memory layout
};


#ifdef USING_R
class PointCloudDataFrame
{
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
  inline size_t kdtree_get_point_count() const { return n_points; }
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const { return coords[dim][idx]; }
  template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }

  // --- Num. points ---
  inline size_t point_count() const { return n_points; }
  inline size_t size()        const { return n_points; }

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

  inline int get_treeid(const size_t idx) const {
    if (!has_treeid()) throw std::runtime_error("Instance segmentation not available in this point cloud");
    return treeid[idx];
  }

  inline void set_treeid(const size_t idx, int v) {
    if (!has_treeid()) throw std::runtime_error("Instance segmentation not available in this point cloud");
    treeid[idx] = v;
  }

  inline double get_pwood(const size_t idx) const {
    if (!has_pwood()) throw std::runtime_error("Wood likelihood data not available in this point cloud");
    return pwood[idx];
  }

  inline void set_pwood(const size_t idx, double v) {
    if (!has_pwood()) throw std::runtime_error("Wood likelihood data not available in this point cloud");
    pwood[idx] = v;
  }

  inline int get_foliage(const size_t idx) const {
    if (!has_foliage()) throw std::runtime_error("Semantic segmentation not available in this point cloud");
    return foliage[idx];
  }

  inline void set_foliage(const size_t idx, int v) {
    if (!has_foliage()) throw std::runtime_error("Semantic segmentation not available in this point cloud");
    foliage[idx] = v;
  }

  inline double get_hag(const size_t idx) const {
    if (!has_hag()) throw std::runtime_error("HAG data not available in this point cloud");
    return hag[idx];
  }

  inline void set_hag(const size_t idx, double v) {
    if (!has_pwood()) throw std::runtime_error("Wood likelihood data not available in this point cloud");
    pwood[idx] = v;
  }

  inline int get_passage(const size_t idx) const {
    if (!passage) throw std::runtime_error("Passage not available in this point cloud");
    return passage[idx];
  }

  inline void set_passage(const size_t idx, int v) {
    if (!passage) throw std::runtime_error("Passage not available in this point cloud");
    passage[idx] = v;
  }

  inline bool is_wood(const size_t idx) const {
    return get_foliage(idx) == 0;
  }

  // --- In-place transforms ---
  void translate(double x, double y, double z) ;
  void scale(double x, double y, double z) ;

  // --- Subset ---
  PointCloudDataFrame subset(const std::vector<bool>& keep, bool xyz_only = false) const;

  // --- Used only in R ----
  explicit PointCloudDataFrame(const Rcpp::DataFrame& df);
  bool linked_to_dataframe() const { return !owns_memory; }

private:
  void cleanup();
  void swap(PointCloudDataFrame& first, PointCloudDataFrame& second) noexcept;
  void init();
  void safe_alloc(size_t n, bool alloc_attrs);

private:
  double* coords[3] = {nullptr, nullptr, nullptr};

  // Optional attributes (nullptr if absent)
  int*    treeid  = nullptr;
  int*    foliage = nullptr;
  int*    passage = nullptr;
  double* hag     = nullptr;
  double* pwood   = nullptr;

  bool owns_memory = false;
  size_t n_points  = 0;
};
#endif

#ifdef USING_R
using PointCloud = PointCloudDataFrame;
#else
using PointCloud = MyPointCloudWrapper;
#endif

#endif

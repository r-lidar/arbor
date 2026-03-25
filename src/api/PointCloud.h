#ifndef ADAPTATOR_H
#define ADAPTATOR_H

#include <vector>
#include <string>
#include <stdexcept>
#include <cstdint>

#ifdef USING_R
#include <Rcpp.h>
#else
#include <span>
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
  inline size_t kdtree_get_point_count() const { return size(); }
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const { return coords[dim][idx]; }
  template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }

  // --- Num. points ---
  inline size_t size()        const { return n_points; }
  inline size_t true_size()   const { return true_n_points; }

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
    if (!has_hag()) throw std::runtime_error("HAG data not available in this point cloud");
    hag[idx] = v;
  }

  inline int get_passage(const size_t idx) const {
    if (!has_passage()) throw std::runtime_error("Passage not available in this point cloud");
    return passage[idx];
  }

  inline void set_passage(const size_t idx, int v) {
    if (!has_passage()) throw std::runtime_error("Passage not available in this point cloud");
    passage[idx] = v;
  }

  inline bool is_wood(const size_t idx) const {
    return get_foliage(idx) == 0;
  }

  inline void set_classification(const size_t idx, int v) {
    if (!has_class()) throw std::runtime_error("Classification data not available in this point cloud");
    classif[idx] = v;
  }

  inline int get_classification(const size_t idx) const {
    if (!has_class()) throw std::runtime_error("Classification data not available in this point cloud");
    return classif[idx];
  }

  inline bool is_ground(const size_t idx) const {
    return get_classification(idx) == 2;
  }

  // -- Partitioning ---
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

  // --- In-place transforms ---
  void translate(double x, double y, double z) ;
  void scale(double x, double y, double z) ;

  // --- Subset ---
  PointCloudDataFrame subset(const std::vector<int>&  indices, bool xyz_only = false) const;
  PointCloudDataFrame subset(const std::vector<bool>& mask,    bool xyz_only = false) const;

  // --- Used only in R ----
  explicit PointCloudDataFrame(const Rcpp::DataFrame& df);
  bool linked_to_dataframe() const { return !owns_memory; }

private:
  void cleanup();
  void swap(PointCloudDataFrame& first, PointCloudDataFrame& second) noexcept;
  void init();
  void safe_alloc(size_t n, bool alloc_attrs);
  void swap_points(size_t i, size_t j)
  {
    for (int d = 0; d < 3; ++d) std::swap(coords[d][i], coords[d][j]);
    if (treeid)  std::swap(treeid[i], treeid[j]);
    if (foliage) std::swap(foliage[i], foliage[j]);
    if (passage) std::swap(passage[i], passage[j]);
    if (classif) std::swap(classif[i], classif[j]);
    if (hag)     std::swap(hag[i], hag[j]);
    if (pwood)   std::swap(pwood[i], pwood[j]);
  }

private:
  double* coords[3] = {nullptr, nullptr, nullptr};

  // Optional attributes (nullptr if absent)
  int*    treeid  = nullptr;
  int*    foliage = nullptr;
  int*    passage = nullptr;
  int*    classif = nullptr;
  double* hag     = nullptr;
  double* pwood   = nullptr;

  bool owns_memory = false;
  size_t n_points  = 0;
  size_t true_n_points = 0;
};
#else
class PointCloudDefault
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
      treeid.push_back(-1);
      foliage.push_back(-1);
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
    q[0] = coords[idx].x;
    q[1] = coords[idx].y;
    q[2] = coords[idx].z;
  }

  inline double get_x(const size_t idx) const { return static_cast<double>(coords[idx].x); }
  inline double get_y(const size_t idx) const { return static_cast<double>(coords[idx].y); }
  inline double get_z(const size_t idx) const { return static_cast<double>(coords[idx].z); }
  inline void   set_x(const size_t idx, double v) { coords[idx].x = v; }
  inline void   set_y(const size_t idx, double v) { coords[idx].y = v; }
  inline void   set_z(const size_t idx, double v) { coords[idx].z = v; }

  // --- Optional attribute access ---
  inline bool has_rgb()     const { return !rgb.empty(); }
  inline bool has_hag()     const { return !hag.empty(); }
  inline bool has_treeid()  const { return !treeid.empty(); }
  inline bool has_pwood()   const { return !pwood.empty(); }
  inline bool has_foliage() const { return !foliage.empty(); }
  inline bool has_passage() const { return !passage.empty(); }
  inline bool has_class()   const { return !classif.empty(); }


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
    if (!has_hag()) throw std::runtime_error("HAG data not available in this point cloud");
    hag[idx] = v;
  }

  inline int get_passage(const size_t idx) const {
    if (!has_passage()) throw std::runtime_error("Passage not available in this point cloud");
    return passage[idx];
  }

  inline void set_passage(const size_t idx, int v) {
    if (!has_passage()) throw std::runtime_error("Passage not available in this point cloud");
    passage[idx] = v;
  }

  inline bool is_wood(const size_t idx) const {
    return get_foliage(idx) == 0;
  }

  inline void set_classification(const size_t idx, int v) {
    if (!has_class()) throw std::runtime_error("Classification data not available in this point cloud");
    classif[idx] = v;
  }

  inline int get_classification(const size_t idx) const {
    if (!has_class()) throw std::runtime_error("Classification data not available in this point cloud");
    return classif[idx];
  }

  inline bool is_ground(const size_t idx) const {
    return get_classification(idx) == 2;
  }

  // -- Partitioning ---
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


  // --- In-place transforms ---
  void translate(double x, double y, double z) ;
  void scale(double x, double y, double z) ;

  // --- Subset ---
  PointCloudDefault subset(const std::vector<int>& indices, bool xyz_only = false) const;
  PointCloudDefault subset(const std::vector<bool>& mask, bool xyz_only = false) const;

  // --- Color ---
  void colorize_trees(bool darken_foliage = false);

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
  std::vector<short> foliage;
  std::vector<short> classif;
  std::vector<int> passage;
  std::vector<float> hag;
  std::vector<float> pwood;

  size_t n_points  = 0;
  size_t true_n_points = 0;
};
#endif

#ifdef USING_R
using PointCloud = PointCloudDataFrame;
#else
using PointCloud = PointCloudDefault;
#endif

#endif

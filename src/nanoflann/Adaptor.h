#ifndef ADAPTATOR_H
#define ADAPTATOR_H

#include <vector>
#include <stdexcept>

#ifdef USING_R
#include <Rcpp.h>
#endif

// ============================================================================
// Virtual base class defining the required interface for all point clouds
// ============================================================================

class PointCloudAdaptorBase
{
public:
  virtual ~PointCloudAdaptorBase() = default;

  // --- Nanoflann KD-tree interface (pure virtual) ---
  virtual size_t kdtree_get_point_count() const = 0;
  virtual double kdtree_get_pt(const size_t idx, const size_t dim) const = 0;
  template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }

  // --- Common interface (pure virtual) ---
  virtual size_t point_count() const = 0;
  virtual size_t size() const = 0;

  // --- Geometry access (pure virtual) ---
  virtual void get_point(const size_t idx, double* q) const = 0;
  virtual double get_x(const size_t idx) const = 0;
  virtual double get_y(const size_t idx) const = 0;
  virtual double get_z(const size_t idx) const = 0;

  // --- Subset (virtual with default implementations) ----
  virtual void subset(const std::vector<bool>& keep)  { throw std::runtime_error("This point cloud has no subset implementation"); };

  // --- In-place transforms (pure virtual) ---
  virtual void translate(double x, double y, double z) = 0;
  virtual void scale(double x, double y, double z) = 0;

  // --- Optional attribute access (virtual with default implementations) ---
  // These can be missing
  virtual bool has_treeid() const { return false; }
  virtual bool has_woodlikelihood() const { return false; }
  virtual bool has_foliage() const { return false; }

  virtual int get_treeid(const size_t idx) const { throw std::runtime_error("Tree ID data not available in this point cloud"); }
  virtual double get_woodlikelihood(const size_t idx) const { throw std::runtime_error("Wood likelihood data not available in this point cloud"); }
  virtual int get_foliage(const size_t idx) const { throw std::runtime_error("Semantic classification not available in this point cloud"); }
  virtual int is_wood(const size_t idx) const { throw std::runtime_error("Semantic classification not available in this point cloud"); }
};

// ============================================================================
// DataFrameAdaptor - Full-featured adaptor with optional attributes
// ============================================================================

#ifdef USING_R
class DataFrameAdaptor : public PointCloudAdaptorBase
{
public:
  // Mandatory coordinates
  const double* coords[3];

  // Optional attributes (nullptr if absent)
  const int*    treeid = nullptr;
  const int*    foliage = nullptr;
  const double* pwood = nullptr;

  bool owns_memory = false;

  size_t n_points;

  DataFrameAdaptor(const Rcpp::DataFrame& df)
  {
    std::vector<std::string> coord_names = {"X", "Y", "Z"};
    std::string treeid_name  = "treeID";
    std::string pwood_name   = "pwood";
    std::string foliage_name = "foliage";

    n_points = df.rows();

    // --- Mandatory coordinates ---
    for (size_t i = 0; i < 3; ++i)
    {
      if (!df.containsElementNamed(coord_names[i].c_str()))
        throw std::runtime_error("Missing mandatory coordinate column: " + coord_names[i]);

      Rcpp::NumericVector col = df[coord_names[i]];
      coords[i] = col.begin();
    }

    // --- Optional attributes ---
    if (df.containsElementNamed(treeid_name.c_str()))
    {
      Rcpp::IntegerVector col = df[treeid_name];
      treeid = col.begin();
    }

    if (df.containsElementNamed(pwood_name.c_str()))
    {
      Rcpp::NumericVector col = df[pwood_name];
      pwood = col.begin();
    }

    if (df.containsElementNamed(foliage_name.c_str()))
    {
      Rcpp::IntegerVector col = df[foliage_name];
      foliage = col.begin();
    }
  }

  ~DataFrameAdaptor() override
  {
    cleanup();
  }

  void cleanup()
  {
    if (owns_memory)
    {
      delete[] coords[0];
      delete[] coords[1];
      delete[] coords[2];
      if (treeid != nullptr) delete[] treeid;
      if (foliage != nullptr) delete[] foliage;
      if (pwood != nullptr) delete[] pwood;

      coords[0] = nullptr;
      coords[1] = nullptr;
      coords[2] = nullptr;
      treeid    = nullptr;
      foliage   = nullptr;
      pwood     = nullptr;
    }
  }

  // --- Nanoflann KD-tree interface ---
  inline size_t kdtree_get_point_count() const override { return n_points; }
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const override
  {
    return coords[dim][idx];
  }

  // --- Num. points ----
  inline size_t point_count() const override { return n_points; }
  inline size_t size() const override { return n_points; }

  // --- Geometry access ---
  inline void get_point(const size_t idx, double* q) const override
  {
    for (size_t d = 0; d < 3; ++d)
      q[d] = coords[d][idx];
  }

  inline double get_x(const size_t idx) const override { return coords[0][idx]; }
  inline double get_y(const size_t idx) const override { return coords[1][idx]; }
  inline double get_z(const size_t idx) const override { return coords[2][idx]; }

  // --- Optional attribute access ---
  inline bool has_treeid() const override { return treeid != nullptr; }
  inline bool has_woodlikelihood() const override { return pwood != nullptr; }
  inline bool has_foliage() const override { return foliage != nullptr; }

  inline int get_treeid(const size_t idx) const override
  {
    if (!treeid) throw std::runtime_error("Tree ID data not available in this point cloud");
    return treeid[idx];
  }

  inline double get_woodlikelihood(const size_t idx) const override
  {
    if (!pwood) throw std::runtime_error("Wood likelihood data not available in this point cloud");
    return pwood[idx];
  }

  inline int get_foliage(const size_t idx) const override
  {
    if (!foliage) throw std::runtime_error("Foliage classification not available in this point cloud");
    return foliage[idx];
  }

  inline int is_wood(const size_t idx) const override
  {
    return get_foliage(idx) == 0;
  }

  // --- In-place transforms ---
  void translate(double x, double y, double z) override
  {
    for (size_t i = 0; i < n_points; ++i)
    {
      if (x != 0) const_cast<double&>(coords[0][i]) -= x;
      if (y != 0) const_cast<double&>(coords[1][i]) -= y;
      if (z != 0) const_cast<double&>(coords[2][i]) -= z;
    }
  }

  void scale(double x, double y, double z) override
  {
    for (size_t i = 0; i < n_points; ++i)
    {
      if (x != 1.0) const_cast<double&>(coords[0][i]) *= x;
      if (y != 1.0) const_cast<double&>(coords[1][i]) *= y;
      if (z != 1.0) const_cast<double&>(coords[2][i]) *= z;
    }
  }

  void subset(const std::vector<bool>& keep) override
  {
    if (keep.size() != n_points)
      throw std::runtime_error("subset mask size mismatch: expected " +  std::to_string(n_points) + " but got " + std::to_string(keep.size()));

    // Count how many points to keep
    size_t new_count = std::count(keep.begin(), keep.end(), true);

    if (new_count == n_points) return;

    // Allocate new storage for coordinates
    double* new_coords[3];
    for (size_t d = 0; d < 3; ++d) new_coords[d] = new double[new_count];

    // Allocate new storage for optional attributes if present
    int* new_treeid = nullptr;
    int* new_foliage = nullptr;
    double* new_pwood = nullptr;

    if (treeid != nullptr)  new_treeid = new int[new_count];
    if (foliage != nullptr) new_foliage = new int[new_count];
    if (pwood != nullptr)   new_pwood = new double[new_count];

    // Copy data for points that are kept
    size_t j = 0;
    for (size_t i = 0; i < n_points; ++i)
    {
      if (keep[i])
      {
        new_coords[0][j] = coords[0][i];
        new_coords[1][j] = coords[1][i];
        new_coords[2][j] = coords[2][i];
        if (treeid != nullptr)  new_treeid[j] = treeid[i];
        if (foliage != nullptr) new_foliage[j] = foliage[i];
        if (pwood != nullptr)   new_pwood[j] = pwood[i];
        j++;
      }
    }

    if (owns_memory) cleanup();
    owns_memory = true;

    // Update pointers to new data
    coords[0] = new_coords[0];
    coords[1] = new_coords[1];
    coords[2] = new_coords[2];
    treeid  = new_treeid;
    pwood   = new_pwood;
    foliage = new_foliage;

    // Update point count
    n_points = new_count;
  }
};
#endif

// ============================================================================
// MatrixAdaptor - Compact matrix-based adaptor
// ============================================================================

#ifdef USING_R
class MatrixAdaptor : public PointCloudAdaptorBase
{
public:
  Rcpp::NumericMatrix& coords;

  MatrixAdaptor(Rcpp::NumericMatrix& m) : coords(m)
  {
    if (coords.ncol() < 3)
      Rcpp::stop("MatrixAdaptor expects at least 3 columns (x, y, z).");
  }

  // --- Nanoflann KD-tree interface ---
  inline size_t kdtree_get_point_count() const override { return coords.nrow(); }
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const override
  {
    return coords(idx, dim);
  }

  // --- Num. points ----
  inline size_t point_count() const override { return coords.nrow(); }
  inline size_t size() const override { return coords.nrow(); }

  // --- Geometry access ---
  inline void get_point(const size_t idx, double* q) const override
  {
    q[0] = coords(idx, 0);
    q[1] = coords(idx, 1);
    q[2] = coords(idx, 2);
  }

  inline double get_x(const size_t idx) const override { return coords(idx, 0); }
  inline double get_y(const size_t idx) const override { return coords(idx, 1); }
  inline double get_z(const size_t idx) const override { return coords(idx, 2); }

  // --- In-place transforms ---
  void translate(double tx, double ty, double tz) override
  {
    if (tx == 0 && ty == 0 && tz == 0) return;
    for (size_t i = 0; i < coords.nrow(); ++i)
    {
      coords(i, 0) -= tx;
      coords(i, 1) -= ty;
      coords(i, 2) -= tz;
    }
  }

  void scale(double sx, double sy, double sz) override
  {
    if (sx == 1.0 && sy == 1.0 && sz == 1.0) return;
    for (size_t i = 0; i < coords.nrow(); ++i)
    {
      coords(i, 0) *= sx;
      coords(i, 1) *= sy;
      coords(i, 2) *= sz;
    }
  }
};
#endif

// ============================================================================
// SimpleAdaptor - Lightweight vector-based adaptor with ID support
// ============================================================================

class SimpleAdaptor : public PointCloudAdaptorBase
{
public:
  struct Point { double x, y, z; int id; };
  std::vector<Point> points;

  // --- Nanoflann KD-tree interface ---
  inline size_t kdtree_get_point_count() const override { return points.size(); }
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const override
  {
    if (dim == 0) return points[idx].x;
    if (dim == 1) return points[idx].y;
    return points[idx].z;
  }

  // --- Num. points ----
  inline size_t point_count() const override { return points.size(); }
  inline size_t size() const override { return points.size(); }

  // --- Geometry access ---
  inline void get_point(const size_t idx, double* q) const override
  {
    q[0] = points[idx].x;
    q[1] = points[idx].y;
    q[2] = points[idx].z;
  }

  inline double get_x(const size_t idx) const override { return points[idx].x; }
  inline double get_y(const size_t idx) const override { return points[idx].y; }
  inline double get_z(const size_t idx) const override { return points[idx].z; }

  // --- In-place transforms ---
  void translate(double tx, double ty, double tz) override
  {
    if (tx == 0 && ty == 0 && tz == 0) return;
    for (auto& pt : points)
    {
      pt.x -= tx;
      pt.y -= ty;
      pt.z -= tz;
    }
  }

  void scale(double sx, double sy, double sz) override
  {
    if (sx == 1.0 && sy == 1.0 && sz == 1.0) return;
    for (auto& pt : points)
    {
      pt.x *= sx;
      pt.y *= sy;
      pt.z *= sz;
    }
  }

  // --- Optional: Access to point ID ---
  inline bool has_treeid() const override { return true; }

  inline int get_treeid(const size_t idx) const override
  {
    return points[idx].id;
  }

  // Optional attributes pwood and foliage use default implementations (not available)
};


using PointCloud = DataFrameAdaptor;

#endif

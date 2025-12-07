#ifndef ADAPTATOR_H
#define ADAPTATOR_H

#include <Rcpp.h>
#include <vector>

class DataFrameAdaptor
{
public:
  std::vector<Rcpp::NumericVector> coords;
  size_t dim;

  DataFrameAdaptor(const Rcpp::DataFrame& df, std::vector<std::string> col_names = {"X", "Y", "Z"})
  {
    dim = col_names.size();
    coords.reserve(dim);
    for (const auto& name : col_names)
      coords.push_back(df[name]);
  }

  inline size_t kdtree_get_point_count() const { return coords[0].size(); }
  inline size_t point_count() const { return coords[0].size(); }
  inline double kdtree_get_pt(const size_t idx, const size_t d) const { return coords[d][idx]; }
  template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
  inline void get_point(const size_t idx, double* q) const
  {
    q[0] = kdtree_get_pt(idx, 0);
    q[1] = kdtree_get_pt(idx, 1);
    q[2] = kdtree_get_pt(idx, 2);
  }
  inline double get_x(const size_t idx) const { return kdtree_get_pt(idx, 0); }
  inline double get_y(const size_t idx) const { return kdtree_get_pt(idx, 1); }
  inline double get_z(const size_t idx) const { return kdtree_get_pt(idx, 2); }
  inline void translate(double x, double y, double z)
  {
    if (x != 0) coords[0] = coords[0] - x;
    if (y != 0) coords[1] = coords[1] - y;
    if (z != 0) coords[2] = coords[2] - z;
  }
  inline void scale(double x, double y, double z)
  {
    if (x != 1.0) coords[0] = coords[0] * x;
    if (y != 1.0) coords[1] = coords[1] * y;
    if (z != 1.0) coords[2] = coords[2] * z;
  }
};

// Compact point cloud adaptor for nanoflann
struct MatrixAdaptor
{
public:
  Rcpp::NumericMatrix& coords;

  explicit MatrixAdaptor(Rcpp::NumericMatrix& m) : coords(m)
  {
    if (coords.ncol() < 3)
      Rcpp::stop("MatrixAdaptor expects at least 3 columns (x, y, z).");
  }

  inline size_t kdtree_get_point_count() const { return coords.nrow(); }
  inline size_t point_count() const { return coords.nrow(); }

  inline double kdtree_get_pt(const size_t idx, const size_t dim) const { return coords(idx, dim); }
  template <class BBOX>  bool kdtree_get_bbox(BBOX&) const { return false; }

  inline void get_point(const size_t idx, double* q) const
  {
    q[0] = coords(idx, 0);
    q[1] = coords(idx, 1);
    q[2] = coords(idx, 2);
  }
  inline double get_x(const size_t idx) { return kdtree_get_pt(idx, 0); }
  inline double get_y(const size_t idx) { return kdtree_get_pt(idx, 1); }
  inline double get_z(const size_t idx) { return kdtree_get_pt(idx, 2); }
  inline void translate(double tx, double ty, double tz)
  {
    if (tx == 0 && ty == 0 && tz == 0) return;
    for (size_t i = 0; i < coords.nrow(); ++i)
    {
      coords(i, 0) -= tx;
      coords(i, 1) -= ty;
      coords(i, 2) -= tz;
    }
  }
  inline void scale(double sx, double sy, double sz)
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

struct SimpleAdaptor
{
  struct Point { double x, y, z; int id; };
  std::vector<Point> points;

  inline size_t kdtree_get_point_count() const { return points.size(); }

  inline double kdtree_get_pt(const size_t idx, const size_t dim) const
  {
    if (dim == 0) return points[idx].x;
    if (dim == 1) return points[idx].y;
    return points[idx].z;
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};


#endif


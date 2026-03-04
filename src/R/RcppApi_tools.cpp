#ifdef USING_R

#include <Rcpp.h>
#include "Grid3D.h"
#include "ransac.h"
#include "arbor.h"

Rcpp::LogicalVector C_homogeneization(Rcpp::DataFrame df, double res, bool hybrid = true)
{
  PointCloud pc(df);
  auto ans = arbor::utils::homogeneization(pc, res, hybrid);
  return(Rcpp::wrap(ans));
}

Rcpp::NumericVector C_anisotropy(Rcpp::DataFrame df,  int k, int ncpu = 1)
{
  PointCloud pc(df);
  auto ans = arbor::utils::anisotropy(pc, k, ncpu);
  return(Rcpp::wrap(ans));
}

Rcpp::IntegerVector C_connected_component(Rcpp::DataFrame df, double res, int connectivity)
{
  PointCloud pc(df);
  Grid3D grid(pc, res);
  return Rcpp::wrap(grid.connected_components(connectivity));
}

Rcpp::LogicalVector C_sor(Rcpp::DataFrame df, unsigned int k, double m, int ncpu = 1)
{
  PointCloud pc(df);
  auto ans = arbor::utils::sor(pc, k, m, ncpu);
  return(Rcpp::wrap(ans));
}


class MatrixAdaptor
{
public:
  Rcpp::NumericMatrix& coords;
  MatrixAdaptor(Rcpp::NumericMatrix& m) : coords(m) { if (coords.ncol() < 3) Rcpp::stop("MatrixAdaptor expects at least 3 columns (x, y, z)."); }
  inline size_t kdtree_get_point_count() const { return coords.nrow(); }
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const { return coords(idx, dim); }
  template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
  inline size_t point_count() const { return coords.nrow(); }
  inline size_t size() const { return coords.nrow(); }
  inline void get_point(const size_t idx, double* q) const { q[0] = coords(idx, 0); q[1] = coords(idx, 1); q[2] = coords(idx, 2); }
  inline double get_x(const size_t idx) const { return coords(idx, 0); }
  inline double get_y(const size_t idx) const { return coords(idx, 1); }
  inline double get_z(const size_t idx) const { return coords(idx, 2); }
};


Rcpp::List ransac_circle_cpp(Rcpp::NumericMatrix x, int num_iterations = 100, double inlier_threshold = 0.01, double early_exit = 1.0)
{
  MatrixAdaptor pc(x);
  RansacCircle rc(num_iterations, inlier_threshold, early_exit);
  for (int i = 0 ; i < pc.point_count() ; i++)
    rc.add_point(pc.get_x(i), pc.get_y(i), pc.get_z(i));
  rc.find_circle();

  std::array<double, 3> center = rc.get_center();
  double radius = rc.get_radius();
  double inlier_pct = rc.get_inlier_percentage();
  double inside_pct = rc.get_inside_percentage();
  double arc_deg = rc.get_arc_coverage();
  const std::vector<int>& inliers = rc.get_inliers();

  // Calculate CFQI (matching R implementation)
  double arc_score = arc_deg / 360.0;
  double score_inside = (inside_pct == 0.0) ? 1.0 : (1.0 - inside_pct);

  double w_arc = 0.4;
  double w_inside = 0.1;
  double w_inlier = 0.4;
  double cfqi = (w_arc * arc_score + w_inlier * inlier_pct + w_inside / score_inside);

  Rcpp::IntegerVector r_inliers = Rcpp::wrap(inliers);

  return Rcpp::List::create(
    Rcpp::Named("center_x") = center[0],
    Rcpp::Named("center_y") = center[1],
    Rcpp::Named("radius") = radius,
    Rcpp::Named("z") = center[2],
    Rcpp::Named("covered_arc_degree") = arc_deg,
    Rcpp::Named("percentage_inlier") = inlier_pct,
    Rcpp::Named("percentage_inside") = inside_pct,
    Rcpp::Named("inliers") = r_inliers+1,
    Rcpp::Named("CFQI") = cfqi
  );
}

#endif

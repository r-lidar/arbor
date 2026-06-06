/**
 * @file RcppApi_tools.cpp
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

#ifdef USING_R

#include <Rcpp.h>
#include "Grid3D.h"
#include "ransac.h"
#include "fitting.h"
#include "allometry.h"
#include "arbor.h"


Rcpp::LogicalVector C_homogeneization(Rcpp::DataFrame df, double res, bool hybrid = true)
{
  PointCloud pc(df);
  auto ans = arbor::utils::homogeneization(pc, res, hybrid);
  return(Rcpp::wrap(ans));
}

Rcpp::NumericVector C_anisotropy(Rcpp::DataFrame df,  int k)
{
  PointCloud pc(df);
  auto ans = arbor::utils::anisotropy(pc, k);
  return(Rcpp::wrap(ans));
}

Rcpp::IntegerVector C_connected_component(Rcpp::DataFrame df, double res, int connectivity)
{
  PointCloud pc(df);
  Grid3D grid(pc, res);
  return Rcpp::wrap(grid.connected_components(connectivity));
}

Rcpp::LogicalVector C_sor(Rcpp::DataFrame df, unsigned int k, double m)
{
  PointCloud pc(df);
  auto ans = arbor::utils::sor(pc, k, m);
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
  for (size_t i = 0 ; i < pc.point_count() ; i++)
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

  Rcpp::IntegerVector r_inliers = Rcpp::wrap(inliers);

  return Rcpp::List::create(
    Rcpp::Named("center_x") = center[0],
    Rcpp::Named("center_y") = center[1],
    Rcpp::Named("radius") = radius,
    Rcpp::Named("z") = center[2],
    Rcpp::Named("covered_arc_degree") = arc_deg,
    Rcpp::Named("percentage_inlier") = inlier_pct,
    Rcpp::Named("percentage_inside") = inside_pct,
    Rcpp::Named("inliers") = r_inliers+1
  );
}

Rcpp::DataFrame allometry(std::string name)
{
  auto model = AllometryDataBase::getAllometry(name);

  std::vector<double> dbh;
  std::vector<double> height;

  for (double h = 0.0; h <= 40.0; h += 0.5)
  {
    double H = model->DBH_vs_H(h);

    height.push_back(h);
    dbh.push_back(H);
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("DBH") = dbh,
    Rcpp::Named("H") = height
  );
}

Rcpp::List fit_circloid_cpp(Rcpp::NumericMatrix x, Rcpp::NumericVector from, Rcpp::NumericVector to, double tolerance, int complexity)
{
  if (x.ncol() != 3) {
    Rcpp::stop("Input matrix must have 3 columns (X, Y, Z)");
  }

  arbor::utils::fitting::Vec3 f = {from[0], from[1], from[2]};
  arbor::utils::fitting::Vec3 t = {to[0], to[1], to[2]};

  arbor::utils::fitting::FittingOrbicular fitter;
  fitter.set_axe(f, t);

  // Add all points
  int n = x.nrow();
  for (int i = 0; i < n; i++) {
    fitter.add_point(x(i, 0), x(i, 1), x(i, 2));
  }

  // Perform fitting
  arbor::utils::fitting::FittingResult result = fitter.fit(tolerance, complexity);

  if (!result.success)
  {
    return Rcpp::List::create(
      Rcpp::Named("success") = false,
      Rcpp::Named("shape_type") = result.shape_type,
      Rcpp::Named("message") = "Fitting failed"
    );
  }

  // Convert inliers to R (1-indexed)
  Rcpp::IntegerVector r_inliers = Rcpp::wrap(result.inlier_indices);
  for (int i = 0; i < r_inliers.size(); i++) { r_inliers[i] += 1; }  // Convert to 1-based indexing

  size_t nnodes = result.nodes.size();
  Rcpp::NumericMatrix nodes(nnodes, 3);

  for (size_t i = 0; i < nnodes; ++i)
  {
    nodes(i,0) = result.nodes[i].x;
    nodes(i,1) = result.nodes[i].y;
    nodes(i,2) = result.nodes[i].z;
  }

  // Build result list based on shape type
  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("success") = result.success,
    Rcpp::Named("radius") = result.radius,
    Rcpp::Named("shape_type") = result.shape_type,
    Rcpp::Named("center_x") = result.center.x,
    Rcpp::Named("center_y") = result.center.y,
    Rcpp::Named("center_z") = result.center.z,
    Rcpp::Named("covered_arc_degree") = result.arc_coverage_deg,
    Rcpp::Named("percentage_inlier") = result.inlier_percentage,
    //Rcpp::Named("percentage_inside") = result.inside_percentge,
    Rcpp::Named("nodes") = nodes,
    Rcpp::Named("inliers") = r_inliers
  );

  return output;
}

namespace arbor::segment{
void fix_small_isolated_low_clusters(PointCloud& las, double res = 0.05, int min_size = 200);
}
void C_fix_small_isolated_low_clusters(Rcpp::DataFrame df, double res = 0.05, int min_size = 200)
{
  PointCloud pc(df);
  arbor::segment::fix_small_isolated_low_clusters(pc, res, min_size);
}

#endif

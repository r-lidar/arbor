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

  // Force data centered on (0,0,0) to normalize geographic coordinates
  if (f.x == 0 && f.y == 0 && t.x == 0 && t.y == 0)
  {
    double mean_x = 0.0;
    double mean_y = 0.0;

    for (int i = 0; i < x.nrow(); ++i)
    {
      mean_x += x(i, 0);
      mean_y += x(i, 1);
    }

    mean_x /= x.nrow();
    mean_y /= x.nrow();

    f.x = mean_x;
    f.y = mean_y;
    t.x = mean_x;
    t.y = mean_y;
  }

  arbor::utils::fitting::CrossSectionFitter fitter;
  fitter.set_axis(f, t);

  // Add all points
  int n = x.nrow();
  for (int i = 0; i < n; i++) {
    fitter.add_point(x(i, 0), x(i, 1), x(i, 2));
  }

  auto strategyFromComplexity = [](int complexity) -> arbor::utils::fitting::FitMode
  {
    switch (complexity)
    {
    case 1:  return arbor::utils::fitting::FitMode::Standard;
    case 2:  return arbor::utils::fitting::FitMode::Standard;
    case 3:  return arbor::utils::fitting::FitMode::Buttress;
    case 4:  return arbor::utils::fitting::FitMode::Full;
    default: return arbor::utils::fitting::FitMode::Basic;
    }
  };

  // Perform fitting
  arbor::utils::fitting::FittingResult result = fitter.fit(tolerance, strategyFromComplexity(complexity));

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

  size_t nnodes = result.contour.size();
  Rcpp::NumericMatrix nodes(nnodes, 3);

  for (size_t i = 0; i < nnodes; ++i)
  {
    nodes(i,0) = result.contour[i].x;
    nodes(i,1) = result.contour[i].y;
    nodes(i,2) = result.contour[i].z;
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

#include "segment_overfitting.h"

namespace arbor::segment {
void fix_small_isolated_low_clusters(PointCloud& las, double res = 0.05, int min_size = 200);
}
void C_fix_small_isolated_low_clusters(Rcpp::DataFrame df, double res = 0.05, int min_size = 200)
{
  PointCloud pc(df);
  arbor::segment::fix_small_isolated_low_clusters(pc, res, min_size);
}

Rcpp::List C_match_instances(Rcpp::DataFrame df)
{
  PointCloud pc(df);
  arbor::segment::MergerConfig config;
  arbor::segment::OverSegmentationResolver merger(config);
  arbor::segment::InstanceMatchResult ans = merger.detect(pc);

  const size_t nc = ans.circles.size();
  Rcpp::NumericVector cx(nc), cy(nc), cz(nc), cr(nc);
  Rcpp::IntegerVector cn(nc);
  for (size_t i = 0; i < nc; ++i)
  {
    cx[i] = ans.circles[i].x;
    cy[i] = ans.circles[i].y;
    cz[i] = ans.circles[i].z;
    cr[i] = ans.circles[i].r;
    cn[i] = static_cast<int>(ans.circles[i].n);
  }
  Rcpp::DataFrame circles = Rcpp::DataFrame::create(
    Rcpp::Named("x") = cx,
    Rcpp::Named("y") = cy,
    Rcpp::Named("z") = cz,
    Rcpp::Named("r") = cr,
    Rcpp::Named("n") = cn,
    Rcpp::Named("stringsAsFactors") = false);

  const size_t nm = ans.matches.size();
  Rcpp::IntegerVector ref(nm), move(nm);
  Rcpp::IntegerVector n(nm);
  Rcpp::NumericVector ref_points(nm), move_points(nm), pct(nm);
  for (size_t i = 0; i < nm; ++i)
  {
    const arbor::segment::MatchRecord& m = ans.matches[i];
    ref[i]         = m.ref;
    move[i]        = m.move;
    n[i]           = static_cast<int>(m.n);
    ref_points[i]  = static_cast<double>(m.ref_points);
    move_points[i] = static_cast<double>(m.move_points);
    pct[i]         = m.weight_ratio();
  }
  Rcpp::DataFrame matching_table = Rcpp::DataFrame::create(
    Rcpp::Named("ref")         = ref,
    Rcpp::Named("move")        = move,
    Rcpp::Named("n")           = n,
    Rcpp::Named("ref_points")  = ref_points,
    Rcpp::Named("move_points") = move_points,
    Rcpp::Named("pct")         = pct,
    Rcpp::Named("stringsAsFactors") = false);

  return Rcpp::List::create(
    Rcpp::Named("circles")        = circles,
    Rcpp::Named("matching_table") = matching_table);
}

void C_merge_instances(Rcpp::List match_result, Rcpp::DataFrame df)
{
  // Extract the 'matching_table' DataFrame from the match_result list
  if (!match_result.containsElementNamed("matching_table"))
  {
    Rcpp::stop("match_result must contain a 'matching_table' element.");
  }
  Rcpp::DataFrame matching_table = Rcpp::as<Rcpp::DataFrame>(match_result["matching_table"]);

  // Extract columns from the matching table
  Rcpp::IntegerVector ref = matching_table["ref"];
  Rcpp::IntegerVector move = matching_table["move"];
  Rcpp::IntegerVector n = matching_table["n"];

  // ref_points/move_points are new columns produced by the current
  // C_match_instances(). Fall back to 0 (i.e. weight_ratio() == 0, which
  // only matters if min_weight_ratio > 0) so a matching_table saved by an
  // older version of this package doesn't hard-fail here.
  const bool has_weights = matching_table.containsElementNamed("ref_points") &&
    matching_table.containsElementNamed("move_points");
  Rcpp::NumericVector ref_points, move_points;
  if (has_weights)
  {
    ref_points  = matching_table["ref_points"];
    move_points = matching_table["move_points"];
  }
  else
  {
    Rcpp::stop("matching_table has no 'ref_points'/'move_points' columns; re-run match_instances() to use min_weight_ratio > 0.");
  }

  // Convert Rcpp vectors into a std::vector<MatchRecord>
  const size_t nm = matching_table.nrows();
  std::vector<arbor::segment::MatchRecord> matches;
  matches.reserve(nm);
  for (size_t i = 0; i < nm; ++i)
  {
    arbor::segment::MatchRecord rec;
    rec.ref  = ref[i];
    rec.move = move[i];
    rec.n    = static_cast<size_t>(n[i]);
    if (has_weights)
    {
      rec.ref_points  = static_cast<size_t>(ref_points[i]);
      rec.move_points = static_cast<size_t>(move_points[i]);
    }
    matches.push_back(rec);
  }

  // Construct PointCloud and set up MergerConfig
  PointCloud pc(df);
  arbor::segment::MergerConfig config;

  // Perform the merge step directly using the extracted matches
  arbor::segment::OverSegmentationResolver merger(config);
  merger.merge(pc, matches);
  return;
}



#endif

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
#include "fitting.h"
#include "allometry.h"
#include "arbor.h"


// For the unit tests of ScopedPartition (tests/testthat/test-partition.R). Partitions df on hag > cut
// and reports what the algorithms see within the scope and after it. Within the scope, writes passage = 1
// for the kept points and passage = 2 for the hidden ones, to check that the attributes written follow the
// points back to their original index. Throws within the scope if fail is true.
Rcpp::List C_scoped_partition(Rcpp::DataFrame df, double cut, bool fail)
{
  PointCloud pc(df);
  const size_t total = pc.size();
  const size_t nqueries = 1000; // knn checks on the first points only: the test clouds are large enough

  // A KD-tree of the whole point cloud before the scope: the partition must invalidate it
  if (total > 0) pc.build_index();

  size_t kept = 0;
  bool kept_above = true;
  bool hidden_below = true;
  bool knn_kept_only = true;
  {
    ScopedPartition<PointCloud> partition(pc, [&](size_t i) { return pc.get_hag(i) > cut; });
    kept = pc.size();

    for (size_t i = 0; i < kept; ++i)
    {
      if (!(pc.get_hag(i) > cut)) kept_above = false;
      pc.set_passage(i, 1);
    }
    for (size_t i = kept; i < partition.total(); ++i)
    {
      if (pc.get_hag(i) > cut) hidden_below = false;
      pc.set_passage(i, 2);
    }

    // The KD-tree is rebuilt on the kept points only: the nearest neighbour of a hidden point is a kept point
    if (kept > 0)
    {
      std::vector<unsigned int> idx;
      std::vector<double> dist;
      double q[3];
      for (size_t i = kept; i < partition.total() && i < kept + nqueries; ++i)
      {
        pc.get_point(i, q);
        pc.knn(q, 1, idx, dist);
        if (idx.empty() || idx[0] >= kept || !(pc.get_hag(idx[0]) > cut)) knn_kept_only = false;
      }
    }

    if (fail) throw std::runtime_error("failure within the scope");
  }

  // The KD-tree rebuilt on the whole point cloud again: every point is its own nearest neighbour (distance 0)
  bool knn_whole = true;
  if (total > 0)
  {
    std::vector<unsigned int> idx;
    std::vector<double> dist;
    double q[3];
    for (size_t i = 0; i < total && i < nqueries; ++i)
    {
      pc.get_point(i, q);
      pc.knn(q, 1, idx, dist);
      if (idx.empty() || dist[0] != 0) knn_whole = false;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("kept") = kept,
    Rcpp::Named("total") = total,
    Rcpp::Named("size_after") = pc.size(),
    Rcpp::Named("kept_above") = kept_above,
    Rcpp::Named("hidden_below") = hidden_below,
    Rcpp::Named("knn_kept_only") = knn_kept_only,
    Rcpp::Named("knn_whole") = knn_whole);
}

Rcpp::LogicalVector C_homogeneization(Rcpp::DataFrame df, double res, bool hybrid = true)
{
  PointCloud pc(df);
  auto ans = arbor::utils::homogeneization(pc, res, hybrid);
  return(Rcpp::wrap(ans));
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

Rcpp::DataFrame allometry_database_cpp()
{
  auto entries = AllometryDataBase::entries();
  R_xlen_t n = entries.size();

  Rcpp::CharacterVector model(n), authors(n), reference(n), country(n);

  for (R_xlen_t i = 0; i < n; i++)
  {
    model[i]     = Rcpp::String(entries[i].name, CE_UTF8);
    authors[i]   = Rcpp::String(entries[i].authors, CE_UTF8);
    reference[i] = Rcpp::String(entries[i].reference, CE_UTF8);
    country[i]   = Rcpp::String(entries[i].country, CE_UTF8);
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("model") = model,
    Rcpp::Named("authors") = authors,
    Rcpp::Named("country") = country,
    Rcpp::Named("reference") = reference,
    Rcpp::Named("stringsAsFactors") = false
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

// Sorts the points of df in place along a Morton (Z-order) curve. Every column follows its point.
void C_spatial_sort(Rcpp::DataFrame df)
{
  PointCloud pc(df);
  pc.spatial_sort();
}

#endif

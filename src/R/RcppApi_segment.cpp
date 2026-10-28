/**
 * @file RcppApi_segment.cpp
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

#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <sstream>

#include "myomp.h"
#include "arbor.h"
#include "GraphBuilder.h"
#include "SeedDetector.h"
#include "RcppApi_wrappers.h"
#include "RcppApi_params.h"


using DF = Rcpp::DataFrame;

void segment_ground_cpp(DF core, Rcpp::List params)
{
  arbor::settings::ArborParameters par = extract_arbor_params(params);
  PointCloud p(core);
  arbor::segment::segment_ground(p, par);
}

void segment_semantic_cpp(DF core, Rcpp::List params)
{
  arbor::settings::ArborParameters par = extract_arbor_params(params);
  PointCloud p(core);
  PointCloud s = arbor::dtm::dtm(p);
  arbor::segment::segment_semantic(p, s, par);
}

void segment_instance_cpp(DF core, DF seeds, Rcpp::List params)
{
  arbor::settings::ArborParameters par = extract_arbor_params(params);
  PointCloud p(core);
  PointCloud s(seeds);
  arbor::segment::segment_instance(p, s, par);
}

DF find_seeds_cpp(DF core, Rcpp::List params)
{
  arbor::settings::ArborParameters par = extract_arbor_params(params);
  PointCloud p(core);
  PointCloud seeds = arbor::seeds::find_seeds(p, par);
  return as_dataframe(seeds);
}

void colorize_trees_cpp(DF core, bool dark)
{
  PointCloud p(core);
  p.colorize_trees(dark);
}

void flag_small_trees_cpp(DF core, double max_height)
{
  PointCloud p(core);
  arbor::segment::flag_small_trees(p, max_height);
}

void flag_buffer_cpp(DF core, Rcpp::Nullable<DF> seeds, double buffer)
{
  PointCloud p(core);
  if (seeds.isNull()) return arbor::segment::flag_buffer(p, buffer);
  PointCloud s(DF(seeds.get()));
  arbor::segment::flag_buffer(p, buffer, &s);
}

void flag_buffer_polygon_cpp(DF core, Rcpp::Nullable<DF> seeds, Rcpp::NumericVector x, Rcpp::NumericVector y)
{
  if (x.size() != y.size()) Rcpp::stop("x and y must have the same length");

  arbor::utils::Polygon::Ring ring(x.size());
  for (R_xlen_t i = 0; i < x.size(); ++i) ring[i] = {x[i], y[i]};
  arbor::utils::Polygon polygon(ring);

  PointCloud p(core);
  if (seeds.isNull()) return arbor::segment::flag_buffer(p, polygon);
  PointCloud s(DF(seeds.get()));
  arbor::segment::flag_buffer(p, polygon, &s);
}

#endif

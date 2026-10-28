/**
 * @file RcppAutoExport.cpp
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

// ! This file exist because Rcpp does not export things that are not in root folder.
// Everything is defined in R/RcppApi_*.cpp
#ifdef USING_R

#include <Rcpp.h>
#include "myomp.h"

//[[Rcpp::export(rng = false)]]
int getThreads()
{
  return omp_get_max_threads();
}

// =======================
// PRE-PROCESSING
// =======================

//[[Rcpp::export(rng = false)]]
Rcpp::List default_arbor_params_cpp();

//[[Rcpp::export(rng = false)]]
Rcpp::LogicalVector C_homogeneization(Rcpp::DataFrame df, double res, bool hybrid = true);

//[[Rcpp::export(rng = false)]]
Rcpp::NumericVector C_anisotropy(Rcpp::DataFrame df,  int k);

// ========================
// SEGMENTATION
// ========================

//[[Rcpp::export(rng = false)]]
void segment_ground_cpp(Rcpp::DataFrame core, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
void segment_semantic_cpp(Rcpp::DataFrame core, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
void segment_instance_cpp(Rcpp::DataFrame core, Rcpp::DataFrame seeds, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
void colorize_trees_cpp(Rcpp::DataFrame core, bool dark);

//[[Rcpp::export(rng = false)]]
void flag_small_trees_cpp(Rcpp::DataFrame core, double max_height);

//[[Rcpp::export(rng = false)]]
void flag_buffer_cpp(Rcpp::DataFrame core, Rcpp::Nullable<Rcpp::DataFrame> seeds, double buffer);

//[[Rcpp::export(rng = false)]]
void flag_buffer_polygon_cpp(Rcpp::DataFrame core, Rcpp::Nullable<Rcpp::DataFrame> seeds, Rcpp::NumericVector x, Rcpp::NumericVector y);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame find_seeds_cpp(Rcpp::DataFrame core, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
Rcpp::List C_match_instances(Rcpp::DataFrame df);

//[[Rcpp::export(rng = false)]]
void C_merge_instances(Rcpp::List match_result, Rcpp::DataFrame df);

// ========================
// QSM
// ========================

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_cpp(Rcpp::DataFrame tree, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
Rcpp::List qsf_cpp(Rcpp::DataFrame scene, double min_height, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
void qsm_write_cpp(Rcpp::DataFrame df, std::string filename, bool binary);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_read_cpp(std::string filename);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_build_and_write_cpp(Rcpp::DataFrame tree, Rcpp::List params, std::string filename);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_dbh_cpp(Rcpp::DataFrame df, double d = 1.3);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_stem_cpp(Rcpp::DataFrame df);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_merchantable_cpp(Rcpp::DataFrame df, double merchantable_radius, double merchantable_length);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_autorepair_cpp(Rcpp::DataFrame df);

//[[Rcpp::export(rng = false)]]
Rcpp::List qsm_mesh_cpp(Rcpp::DataFrame df, int resolution);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame allometry(std::string name);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame allometry_database_cpp();

// ========================
// QSF
// ========================

//[[Rcpp::export(rng = false)]]
void qsf_write_cpp(Rcpp::List x, std::string dir, std::string format, bool binary);

//[[Rcpp::export(rng = false)]]
Rcpp::List qsf_read_cpp(std::string filename);

// ========================
// FITTING
// ========================

//[[Rcpp::export(rng = false)]]
Rcpp::List fit_circloid_cpp(Rcpp::NumericMatrix x, Rcpp::NumericVector from = Rcpp::NumericVector::create(0,0,0), Rcpp::NumericVector to = Rcpp::NumericVector::create(0,0,1), double tolerance = 0.02, int complexity = 1);

// ========================
// INVENTORY MERGING
// ========================

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector aim_disentangle_ids_cpp(Rcpp::DataFrame las, Rcpp::DataFrame inventory);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame aim_matching_table_cpp(Rcpp::DataFrame las, Rcpp::DataFrame inventory);

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector aim_reassign_ids_cpp(Rcpp::DataFrame las, Rcpp::DataFrame matching_table);

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector aim_resolve_multi_matching_cpp(Rcpp::DataFrame las, Rcpp::DataFrame matching_table, Rcpp::DataFrame inventory);

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector aim_merge_inventory_cpp(Rcpp::DataFrame las, Rcpp::DataFrame inventory);

//[[Rcpp::export(rng = false)]]
Rcpp::List aim_fit_tree_cpp(Rcpp::DataFrame las, Rcpp::DataFrame tree);

// ========================
// EXPERIMENTAL
// ========================

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_distances_cpp(Rcpp::DataFrame qsm_df, Rcpp::DataFrame pts_df);

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector extract_tree_context_cpp(Rcpp::DataFrame las, int tree_id,  bool exclude_tree = false, int k = 10);

#endif

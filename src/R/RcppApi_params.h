/**
 * @file RcppApi_params.h
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

#ifndef APIPARAM_H
#define APIPARAM_H

#ifdef USING_R

#include "arbor.h"

inline void assert_exists(const Rcpp::List& p, const char* name)
{
  if (!p.containsElementNamed(name))
  {
    Rcpp::stop("Invalid parameters: missing '%s'", name);
  }
}

// Helper function to extract parameters from R list
arbor::settings::GlobalParameters extract_global_params(const Rcpp::List& params);
arbor::settings::GraphParameters extract_pathfinder_params(Rcpp::List params);
arbor::settings::SemanticParameters extract_semantic_params(const Rcpp::List& params);
arbor::settings::WoodlikelihoodParameters extract_likelihood_params(const Rcpp::List& params);
arbor::settings::SeedParameters extract_seeds_params(const Rcpp::List& params);
arbor::settings::QsmParameters extract_qsm_params(const Rcpp::List& params);
arbor::settings::ArborParameters extract_arbor_params(const Rcpp::List& params);

// Helper function to convert to R list
Rcpp::List global_to_list(const arbor::settings::GlobalParameters& s);
Rcpp::List likelihood_to_list(const arbor::settings::WoodlikelihoodParameters& s);
Rcpp::List semantic_to_list(const arbor::settings::SemanticParameters& s);
Rcpp::List seeds_to_list(const arbor::settings::SeedParameters& s);
Rcpp::List qsm_to_list(const arbor::settings::QsmParameters& s);
Rcpp::List graph_to_pathfinder_list(const arbor::settings::GraphParameters& g);
Rcpp::List graph_to_instance_list(const arbor::settings::GraphParameters& g);
Rcpp::List default_arbor_params_cpp(const arbor::settings::ArborParameters& params);

void Rlogger(const std::string& msg);

#endif

#endif

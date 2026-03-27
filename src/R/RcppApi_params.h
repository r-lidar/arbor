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

void Rlogger(const std::string& msg);

#endif

#endif

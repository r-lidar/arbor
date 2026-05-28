/**
 * @file RcppApi_params.cpp
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

#include "RcppApi_params.h"
#include "services.h"

#include <chrono>
#include <iomanip>
#include <string>
#include <sstream>

// ---- Logger -----------------------------------------------------------

void Rlogger(const std::string& msg)
{
  auto now = std::chrono::system_clock::now();
  auto time_t_now = std::chrono::system_clock::to_time_t(now);
  std::tm* tm_now = std::localtime(&time_t_now);

  std::ostringstream oss;
  oss << "["
      << std::setfill('0') << std::setw(2) << tm_now->tm_hour << ":"
      << std::setfill('0') << std::setw(2) << tm_now->tm_min  << ":"
      << std::setfill('0') << std::setw(2) << tm_now->tm_sec
      << "] " << msg;

  Rcpp::Rcout << oss.str() << std::endl;
}

namespace {
struct RServicesInit {
  RServicesInit() { ServiceLocator::register_logger(Rlogger); }
};
  static RServicesInit r_services_init;
}

// ---- Visitors — written once, never again ----------------------------

struct RExtractor
{
  const Rcpp::List& p;

  template<typename T>
  void operator()(const char* name, T& field)
  {
    assert_exists(p, name);
    field = Rcpp::as<T>(p[name]);
  }
};

struct RSerializer
{
  Rcpp::List out;

  template<typename T>
  void operator()(const char* name, T& field) { out[name] = field; }
};

// ---- Generic helpers — written once, never again ---------------------

template<typename S>
S extract(const Rcpp::List& params, const char* key)
{
  assert_exists(params, key);
  Rcpp::List p = params[key];
  S s;
  RExtractor v{p};
  s.visit(v);
  return s;
}

template<typename S>
Rcpp::List to_list(S& s)
{
  RSerializer v;
  s.visit(v);
  return v.out;
}

// ---- Extract ----------------------------------------------------------

arbor::settings::GlobalParameters extract_global_params(const Rcpp::List& params)
{
  return extract<arbor::settings::GlobalParameters>(params, "global");
}

arbor::settings::WoodlikelihoodParameters extract_likelihood_params(const Rcpp::List& params)
{
  return extract<arbor::settings::WoodlikelihoodParameters>(params, "woodlikelihood");
}

arbor::settings::SemanticParameters extract_semantic_params(const Rcpp::List& params)
{
  return extract<arbor::settings::SemanticParameters>(params, "semantic");
}

arbor::settings::QsmParameters extract_qsm_params(const Rcpp::List& params)
{
  return extract<arbor::settings::QsmParameters>(params, "qsm");
}

arbor::settings::SeedParameters extract_seeds_params(const Rcpp::List& params)
{
  assert_exists(params, "seed");
  Rcpp::List p = params["seed"];
  arbor::settings::SeedParameters s;
  RExtractor v{p};
  s.visit(v);                                                    // handles scalar fields
  assert_exists(p, "slice_at");
  s.slice_at = Rcpp::as<std::vector<double>>(p["slice_at"]);    // special: vector
  return s;
}

arbor::settings::GraphParameters extract_pathfinder_params(Rcpp::List params)
{
  assert_exists(params, "path_finder");
  Rcpp::List p = params["path_finder"];
  arbor::settings::GraphParameters s;
  RExtractor pf{p};
  s.visit_pathfinder(pf);                                        // k, decimation, space_res ...
  assert_exists(p, "penalty");
  s.angle_penalty = Rcpp::as<std::vector<float>>(p["penalty"]); // special: computed default

  assert_exists(params, "instance");
  p = params["instance"];
  RExtractor inst{p};
  s.visit_instance(inst);                                        // wood2leaf, leaf2leaf, wood2wood
  return s;
}

arbor::settings::ArborParameters extract_arbor_params(const Rcpp::List& params)
{
  arbor::settings::ArborParameters s;
  s.global         = extract_global_params(params);
  s.woodlikelihood = extract_likelihood_params(params);
  s.pathfinder     = extract_pathfinder_params(params);
  s.semantic       = extract_semantic_params(params);
  s.seeds          = extract_seeds_params(params);
  s.qsm            = extract_qsm_params(params);
  return s;
}

// ---- Serialize --------------------------------------------------------

Rcpp::List global_to_list(arbor::settings::GlobalParameters& s)
{ return to_list(s); }

Rcpp::List likelihood_to_list(arbor::settings::WoodlikelihoodParameters& s)
{ return to_list(s); }

Rcpp::List semantic_to_list(arbor::settings::SemanticParameters& s)
{ return to_list(s); }

Rcpp::List qsm_to_list(arbor::settings::QsmParameters& s)
{ return to_list(s); }

Rcpp::List seeds_to_list(arbor::settings::SeedParameters& s)
{
  RSerializer v;
  s.visit(v);                                   // handles scalar fields
  v.out["slice_at"] = Rcpp::wrap(s.slice_at);  // special: needs wrap
  return v.out;
}

Rcpp::List graph_to_pathfinder_list(arbor::settings::GraphParameters& s)
{
  RSerializer v;
  s.visit_pathfinder(v);
  v.out["penalty"] = Rcpp::wrap(s.angle_penalty); // special: needs wrap
  return v.out;
}

Rcpp::List graph_to_instance_list(arbor::settings::GraphParameters& s)
{
  RSerializer v;
  s.visit_instance(v);
  return v.out;
}

Rcpp::List default_arbor_params_cpp()
{
  arbor::settings::ArborParameters p;

  return Rcpp::List::create(
    Rcpp::Named("global")         = global_to_list(p.global),
    Rcpp::Named("woodlikelihood") = likelihood_to_list(p.woodlikelihood),
    Rcpp::Named("path_finder")    = graph_to_pathfinder_list(p.pathfinder),
    Rcpp::Named("instance")       = graph_to_instance_list(p.pathfinder),
    Rcpp::Named("semantic")       = semantic_to_list(p.semantic),
    Rcpp::Named("seed")           = seeds_to_list(p.seeds),
    Rcpp::Named("qsm")            = qsm_to_list(p.qsm)
  );
}

#endif

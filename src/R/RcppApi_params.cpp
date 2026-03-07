#ifdef USING_R

#include "RcppApi_params.h"

#include <chrono>
#include <string>
#include <sstream>

void Rlogger(const std::string& msg)
{
  auto now = std::chrono::system_clock::now();
  auto time_t_now = std::chrono::system_clock::to_time_t(now);
  std::tm* tm_now = std::localtime(&time_t_now);

  std::ostringstream oss;
  oss << "["
      << std::setfill('0') << std::setw(2) << tm_now->tm_hour << ":"
      << std::setfill('0') << std::setw(2) << tm_now->tm_min << ":"
      << std::setfill('0') << std::setw(2) << tm_now->tm_sec
      << "] " << msg;

  Rcpp::Rcout << oss.str() << std::endl;
}

// Helper function to extract parameters from R list
arbor::settings::GraphParameters extract_pathfinder_params(Rcpp::List params)
{
  assert_exists(params, "path_finder");
  Rcpp::List p = params["path_finder"];
  assert_exists(p, "k_neighborhood_connectivity");
  assert_exists(p, "k_seed_connectivity");
  assert_exists(p, "decimation");
  assert_exists(p, "space_res");
  assert_exists(p, "max_gap");
  assert_exists(p, "penalty");
  assert_exists(p, "distance_power");

  arbor::settings::GraphParameters g;
  g.k = Rcpp::as<int>(p["k_neighborhood_connectivity"]);
  g.k_seed = Rcpp::as<int>(p["k_seed_connectivity"]);
  g.decimation = Rcpp::as<double>(p["decimation"]);
  g.space_res = Rcpp::as<double>(p["space_res"]);
  g.max_gap = Rcpp::as<double>(p["max_gap"]);
  g.power = Rcpp::as<double>(p["distance_power"]);
  g.angle_penalty = Rcpp::as<std::vector<float>>(p["penalty"]);


  assert_exists(params, "instance");
  p = params["instance"];
  assert_exists(p, "wood2leaf_factor");
  assert_exists(p, "leaf2leaf_factor");
  assert_exists(p, "wood2wood_factor");

  g.wood2leaf = Rcpp::as<double>(p["wood2leaf_factor"]);
  g.leaf2leaf = Rcpp::as<double>(p["leaf2leaf_factor"]);
  g.wood2wood = Rcpp::as<double>(p["wood2wood_factor"]);

  return g;
}

arbor::settings::SemanticParameters extract_semantic_params(const Rcpp::List& params)
{
  assert_exists(params, "semantic");
  Rcpp::List p = params["semantic"];

  assert_exists(p, "min_passage");
  assert_exists(p, "high_pwood_threshold");
  assert_exists(p, "medium_pwood_thresold");
  assert_exists(p, "connected_components_res");
  assert_exists(p, "connected_components_min");
  assert_exists(p, "wood_assignation_k");
  assert_exists(p, "wood_assignation_dist");
  assert_exists(p, "wood_extra_reasignation_k");
  assert_exists(p, "wood_extra_reasignation_dist");
  assert_exists(p, "medium_pwood_sor_k");
  assert_exists(p, "medium_pwood_sor_m");
  assert_exists(p, "ground_res");

  arbor::settings::SemanticParameters s;
  s.min_passage = Rcpp::as<int>(p["min_passage"]);
  s.high_pwood_threshold   = Rcpp::as<double>(p["high_pwood_threshold"]);
  s.medium_pwood_threshold = Rcpp::as<double>(p["medium_pwood_thresold"]);
  s.connected_components_res = Rcpp::as<double>(p["connected_components_res"]);
  s.connected_components_min = Rcpp::as<int>(p["connected_components_min"]);
  s.wood_assignation_k = Rcpp::as<int>(p["wood_assignation_k"]);
  s.wood_assignation_dist = Rcpp::as<double>(p["wood_assignation_dist"]);
  s.wood_extra_reasignation_k = Rcpp::as<int>(p["wood_extra_reasignation_k"]);
  s.wood_extra_reasignation_dist = Rcpp::as<double>(p["wood_extra_reasignation_dist"]);
  s.medium_pwood_sor_k = Rcpp::as<int>(p["medium_pwood_sor_k"]);
  s.medium_pwood_sor_m = Rcpp::as<double>(p["medium_pwood_sor_m"]);
  s.ground_res = Rcpp::as<double>(p["ground_res"]);

  return s;
}

arbor::settings::WoodlikelihoodParameters extract_likelihood_params(const Rcpp::List& params)
{
  assert_exists(params, "woodlikelihood");
  Rcpp::List p = params["woodlikelihood"];

  assert_exists(p, "k");

  arbor::settings::WoodlikelihoodParameters s;
  return s;
}


arbor::settings::SeedParameters extract_seeds_params(const Rcpp::List& params)
{
  assert_exists(params, "seed");
  Rcpp::List p = params["seed"];

  assert_exists(p, "slice_at");
  assert_exists(p, "slice_thickness");
  assert_exists(p, "min_passage");
  assert_exists(p, "safe_zone");

  arbor::settings::SeedParameters s;
  s.min_passage = Rcpp::as<int>(p["min_passage"]);
  s.slice_thickness   = Rcpp::as<double>(p["slice_thickness"]);
  s.slice_at = Rcpp::as<std::vector<double>>(p["slice_at"]);
  s.safe_zone = Rcpp::as<double>(p["safe_zone"]);

  return s;
}

arbor::settings::QsmParameters extract_qsm_params(const Rcpp::List& params)
{
  assert_exists(params, "qsm");
  Rcpp::List p = params["qsm"];

  assert_exists(p, "step");
  assert_exists(p, "cl_dist");
  assert_exists(p, "max_d");
  assert_exists(p, "apex");
  assert_exists(p, "smooth_steps");
  assert_exists(p, "smooth_lambda");
  assert_exists(p, "smooth_mu");


  arbor::settings::QsmParameters s;
  s.step   = Rcpp::as<double>(p["step"]);
  s.cl_dist = Rcpp::as<double>(p["cl_dist"]);
  s.max_d = Rcpp::as<double>(p["max_d"]);
  s.apex = Rcpp::as<double>(p["apex"]);
  s.smooth_steps = Rcpp::as<int>(p["smooth_steps"]);
  s.smooth_lambda = Rcpp::as<double>(p["smooth_lambda"]);
  s.smooth_mu = Rcpp::as<double>(p["smooth_mu"]);

  return s;
}

arbor::settings::ArborParameters extract_arbor_params(const Rcpp::List& params)
{
  arbor::settings::GraphParameters gp = extract_pathfinder_params(params);
  arbor::settings::SemanticParameters sp = extract_semantic_params(params);
  arbor::settings::WoodlikelihoodParameters wp = extract_likelihood_params(params);
  arbor::settings::SeedParameters ep = extract_seeds_params(params);
  arbor::settings::QsmParameters qp = extract_qsm_params(params);

  arbor::settings::ArborParameters s;
  s.pathfinder = gp;
  s.semantic = sp;
  s.woodlikelihood = wp;
  s.seeds = ep;
  s.qsm = qp;

  return s;
}

#endif

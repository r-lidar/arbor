#include <Rcpp.h>

#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <sstream>

#include "myomp.h"
#include "api.h"

// Type aliases for clarity
using GraphPtr = Rcpp::XPtr<Graph>;
using GraphCache = std::pair<DistanceVector, PredecessorMap>;
using DF = Rcpp::DataFrame;

auto logger()
{
  return [](const std::string& msg)
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
  };
}


// Inline helper to check required list elements
inline void assert_exists(const Rcpp::List& p, const char* name)
{
  if (!p.containsElementNamed(name))
  {
    Rcpp::stop("Invalid parameters: missing '%s'", name);
  }
}

// Helper function to extract parameters from R list
GraphBuilderParams extract_pathfinder_params(Rcpp::List params)
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

  GraphBuilderParams g;
  g.k = Rcpp::as<int>(p["k_neighborhood_connectivity"]);
  g.k_seed = Rcpp::as<int>(p["k_seed_connectivity"]);
  g.decimation = Rcpp::as<double>(p["decimation"]);
  g.space_res = Rcpp::as<double>(p["space_res"]);
  g.max_gap = Rcpp::as<double>(p["max_gap"]);
  g.power = Rcpp::as<double>(p["distance_power"]);
  g.angle_penalty = Rcpp::as<std::vector<float>>(p["penalty"]);

  return g;
}

SemanticParams extract_semantic_params(const Rcpp::List& params)
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

  SemanticParams s;
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


Rcpp::IntegerVector segment_instance_cpp(DF core, DF seeds, Rcpp::List params)
{
  GraphBuilderParams gparams = extract_pathfinder_params(params);
  PointCloud p(core);
  PointCloud s(seeds);
  std::vector<int> ans = segment_instance(p, s, gparams, logger());
  for (auto& id : ans) { if (id == -1) id = NA_INTEGER; }
  return Rcpp::IntegerVector(ans.begin(), ans.end());
}

Rcpp::IntegerVector accumulate_passages_cpp(DF core, DF gnd, Rcpp::List params)
{
  GraphBuilderParams gparams = extract_pathfinder_params(params);
  PointCloud p(core);
  PointCloud s(gnd);
  std::vector<int> ans = accumulate_passages(p, s, gparams, logger());
  return Rcpp::IntegerVector(ans.begin(), ans.end());
}

Rcpp::LogicalVector assign_wood_from_passage_cpp(DF core, Rcpp::List params)
{
  SemanticParams sparams = extract_semantic_params(params);
  PointCloud p(core);
  std::vector<bool> ans = assign_wood_from_passage(p, sparams, logger());
  return Rcpp::LogicalVector(ans.begin(), ans.end());
}

Rcpp::LogicalVector assign_wood_from_high_likelihood_cpp(DF core, Rcpp::List params)
{
  SemanticParams sparams = extract_semantic_params(params);
  PointCloud p(core);
  std::vector<bool> ans = assign_wood_from_high_likelihood(p, sparams, logger());
  return Rcpp::LogicalVector(ans.begin(), ans.end());
}

Rcpp::LogicalVector assign_wood_from_medium_likelihood_cpp(DF core, Rcpp::List params)
{
  SemanticParams sparams = extract_semantic_params(params);
  PointCloud p(core);
  std::vector<bool> ans = assign_wood_from_medium_likelihood(p, sparams, logger());
  return Rcpp::LogicalVector(ans.begin(), ans.end());
}

Rcpp::LogicalVector assign_wood_from_wood_dilatation_cpp(DF core, Rcpp::List params)
{
  SemanticParams sparams = extract_semantic_params(params);
  PointCloud p(core);
  std::vector<bool> ans = assign_wood_from_wood_dilatation(p, sparams, logger());
  return Rcpp::LogicalVector(ans.begin(), ans.end());
}

SEXP build_semantic_graph(DF dec, DF targets, DF gnd, Rcpp::List params)
{
  GraphBuilderParams gparams = extract_pathfinder_params(params);

  PointCloud core(dec);
  PointCloud trgt(targets);
  PointCloud ground(gnd);

  Graph* graph = build_semantic_graph(core, trgt, ground, gparams);
  GraphPtr ptr(graph, true);
  return ptr;
}

SEXP build_instance_graph(DF dec, DF seed, Rcpp::List params)
{
  GraphBuilderParams gparams = extract_pathfinder_params(params);

  PointCloud core(dec);
  PointCloud seeds(seed);

  Graph* graph = build_instance_graph(core, seeds, gparams);

  GraphPtr ptr(graph, true);
  return ptr;
}

Rcpp::IntegerVector accumulate_passages_old(SEXP graph_ptr, int start_node, Rcpp::IntegerVector goal_nodes, int num_points)
{
  GraphPtr graph(graph_ptr);
  const int n_goals = goal_nodes.size();

  // Global count vector
  std::vector<int> passage(num_points, 0);

  // Precompute distances for fast access
  GraphCache cache = graph->compute_distances(start_node);

  // Parallel loop over goal nodes
  #pragma omp parallel
  {
    std::vector<int> local_passage(num_points, 0);  // thread-local counts

    #pragma omp for schedule(dynamic, 100)
    for (int i = 0; i < n_goals; ++i)
    {
      NodeId goal  = goal_nodes[i];

      auto [path, cost] = graph->findPath(start_node, goal, cache);

      for (size_t j = 0; j < path.size(); ++j)
      {
        NodeId id = path[j];
        if (id >= 0 && id < num_points)
          local_passage[id] += 1;
      }
    }

    // Merge results into global passage safely
    #pragma omp critical
    {
      for (int i = 0; i < num_points; ++i)
        passage[i] += local_passage[i];
    }
  }

return Rcpp::wrap(passage);
}

/*
 * For a set query nodes (usually seed nodes) computes for each node of the graph
 * the distance to the closest query node and return the id of the
 * closest nodes
 */
Rcpp::IntegerVector find_closest_node(SEXP graph_ptr, Rcpp::IntegerVector ids)
{
  GraphPtr graph(graph_ptr);

  // Convert R vector of ground nodes to std::vector<NodeId>
  NodeIDs node_ids(ids.begin(), ids.end());

  std::vector<double> distances;
  NodeIDs closest_nodeids;

  // Run the optimized multi-source Dijkstra
  graph->shortest_paths_from_node(node_ids, distances, closest_nodeids);

  // Return as R list
  return Rcpp::IntegerVector(closest_nodeids.begin(), closest_nodeids.end());
}

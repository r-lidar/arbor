#include <Rcpp.h>

#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <sstream>

#include "myomp.h"
#include "api.h"
#include "GraphBuilder.h"
#include "SeedDetector.h"
#include "Rwrappers.h"


// Type aliases for clarity
using GraphPtr = Rcpp::XPtr<Graph>;
using GraphCache = std::pair<Graph::DistanceVector, Graph::PredecessorMap>;
using DF = Rcpp::DataFrame;

std::vector<int>  accumulate_passages(const PointCloud& core, const PointCloud& ground, const GraphParameters& params, const Logger& logger = [](const std::string&) {});
std::vector<bool> assign_wood_from_passage(const PointCloud& pc, const SemanticParameters& params, const Logger& logger = [](const std::string&) {});
std::vector<bool> assign_wood_from_high_likelihood(const PointCloud& pc, const SemanticParameters& params, const Logger& logger = [](const std::string&) {});
std::vector<bool> assign_wood_from_medium_likelihood(const PointCloud& pc, const SemanticParameters& params, const Logger& logger = [](const std::string&) {});
std::vector<bool> assign_wood_from_wood_dilatation(const PointCloud& pc, const SemanticParameters& params, const Logger& logger = [](const std::string&) {});
Graph* build_semantic_graph(const PointCloud& core, const PointCloud& target, const PointCloud& gnd, const GraphParameters& params);
Graph* build_instance_graph(const PointCloud& core, const PointCloud& seeds, const GraphParameters& params);

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
GraphParameters extract_pathfinder_params(Rcpp::List params)
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

  GraphParameters g;
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

SemanticParameters extract_semantic_params(const Rcpp::List& params)
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

  SemanticParameters s;
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

WoodlikelihoodParameters extract_likelihood_params(const Rcpp::List& params)
{
  assert_exists(params, "woodlikelihood");
  Rcpp::List p = params["woodlikelihood"];

  assert_exists(p, "k");

  WoodlikelihoodParameters s;
  return s;
}


SeedParameters extract_seeds_params(const Rcpp::List& params)
{
  assert_exists(params, "seed");
  Rcpp::List p = params["seed"];

  assert_exists(p, "slice_at");
  assert_exists(p, "slice_thickness");
  assert_exists(p, "min_passage");
  assert_exists(p, "safe_zone");

  SeedParameters s;
  s.min_passage = Rcpp::as<int>(p["min_passage"]);
  s.slice_thickness   = Rcpp::as<double>(p["slice_thickness"]);
  s.slice_at = Rcpp::as<std::vector<double>>(p["slice_at"]);
  s.safe_zone = Rcpp::as<double>(p["safe_zone"]);

  return s;
}

ArborParameters extract_arbor_params(const Rcpp::List& params)
{
  GraphParameters gp = extract_pathfinder_params(params);
  SemanticParameters sp = extract_semantic_params(params);
  WoodlikelihoodParameters wp = extract_likelihood_params(params);
  SeedParameters ep = extract_seeds_params(params);

  ArborParameters s;
  s.pathfinder = gp;
  s.semantic = sp;
  s.woodlikelihood = wp;
  s.seeds = ep;
  return s;
}



void segment_semantic_cpp(DF core, DF ground, Rcpp::List params)
{
  ArborParameters par = extract_arbor_params(params);
  PointCloud p(core);
  PointCloud s(ground);
  segment_semantic(p, s, par, logger());
}

void segment_instance_cpp(DF core, DF seeds, Rcpp::List params)
{
  ArborParameters par = extract_arbor_params(params);
  PointCloud p(core);
  PointCloud s(seeds);
  segment_instance(p, s, par, logger());
  for (size_t i = 0 ; i < p.size() ; i++) {
    if (p.get_treeid(i) == -1) {
      p.set_treeid(i, NA_INTEGER);
    }
  }
}

DF find_seeds_cpp(DF core, Rcpp::List params)
{
  ArborParameters par = extract_arbor_params(params);
  PointCloud p(core);
  PointCloud seeds = find_seeds(p, par, logger());
  return as_dataframe(seeds);
}

Rcpp::IntegerVector accumulate_passages_cpp(DF core, DF gnd, Rcpp::List params)
{
  GraphParameters gparams = extract_pathfinder_params(params);
  PointCloud p(core);
  PointCloud s(gnd);
  std::vector<int> ans = accumulate_passages(p, s, gparams, logger());
  return Rcpp::IntegerVector(ans.begin(), ans.end());
}

Rcpp::LogicalVector assign_wood_from_passage_cpp(DF core, Rcpp::List params)
{
  SemanticParameters sparams = extract_semantic_params(params);
  PointCloud p(core);
  std::vector<bool> ans = assign_wood_from_passage(p, sparams, logger());
  return Rcpp::LogicalVector(ans.begin(), ans.end());
}

Rcpp::LogicalVector assign_wood_from_high_likelihood_cpp(DF core, Rcpp::List params)
{
  SemanticParameters sparams = extract_semantic_params(params);
  PointCloud p(core);
  std::vector<bool> ans = assign_wood_from_high_likelihood(p, sparams, logger());
  return Rcpp::LogicalVector(ans.begin(), ans.end());
}

Rcpp::LogicalVector assign_wood_from_medium_likelihood_cpp(DF core, Rcpp::List params)
{
  SemanticParameters sparams = extract_semantic_params(params);
  PointCloud p(core);
  std::vector<bool> ans = assign_wood_from_medium_likelihood(p, sparams, logger());
  return Rcpp::LogicalVector(ans.begin(), ans.end());
}

Rcpp::LogicalVector assign_wood_from_wood_dilatation_cpp(DF core, Rcpp::List params)
{
  SemanticParameters sparams = extract_semantic_params(params);
  PointCloud p(core);
  std::vector<bool> ans = assign_wood_from_wood_dilatation(p, sparams, logger());
  return Rcpp::LogicalVector(ans.begin(), ans.end());
}

SEXP build_semantic_graph(DF dec, DF targets, DF gnd, Rcpp::List params)
{
  GraphParameters p = extract_pathfinder_params(params);

  PointCloud core(dec);
  PointCloud trgt(targets);
  PointCloud ground(gnd);

  Graph* graph = build_semantic_graph(core, trgt, ground, p);
  GraphPtr ptr(graph, true);
  return ptr;
}

SEXP build_instance_graph(DF dec, DF seed, Rcpp::List params)
{
  GraphParameters p = extract_pathfinder_params(params);

  PointCloud core(dec);
  PointCloud seeds(seed);

  Graph* graph = build_instance_graph(core, seeds, p);

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
  Graph::GraphCache cache = graph->compute_distances(start_node);

  // Parallel loop over goal nodes
  #pragma omp parallel
  {
    std::vector<int> local_passage(num_points, 0);  // thread-local counts

    #pragma omp for schedule(dynamic, 100)
    for (int i = 0; i < n_goals; ++i)
    {
      Graph::NodeId goal  = goal_nodes[i];

      auto [path, cost] = graph->findPath(start_node, goal, cache);

      for (size_t j = 0; j < path.size(); ++j)
      {
        Graph::NodeId id = path[j];
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
  Graph::NodeIDs node_ids(ids.begin(), ids.end());

  std::vector<double> distances;
  Graph::NodeIDs closest_nodeids;

  // Run the optimized multi-source Dijkstra
  graph->shortest_paths_from_node(node_ids, distances, closest_nodeids);

  // Return as R list
  return Rcpp::IntegerVector(closest_nodeids.begin(), closest_nodeids.end());
}


Rcpp::DataFrame generate_cage_cpp(Rcpp::DataFrame circles, double decimation)
{
  // Convert input: R DataFrame -> C++ vector of Circles
  Rcpp::NumericVector X = circles["X"];
  Rcpp::NumericVector Y = circles["Y"];
  Rcpp::NumericVector Z = circles["Z"];
  Rcpp::NumericVector R = circles["R"];
  Rcpp::IntegerVector id = circles["id"];

  int n = X.size();

  if (n == 0) {
    return Rcpp::DataFrame::create(
      Rcpp::Named("X") = Rcpp::NumericVector(),
      Rcpp::Named("Y") = Rcpp::NumericVector(),
      Rcpp::Named("Z") = Rcpp::NumericVector()
    );
  }

  std::vector<SeedDetectorGeometries::Circle> circle_vec;
  circle_vec.reserve(n);
  for (int i = 0; i < n; ++i) {

    circle_vec.push_back(SeedDetectorGeometries::Circle(X[i], Y[i], Z[i], R[i], id[i]));
  }

  // Call pure C++ function
  std::vector<SeedDetectorGeometries::Point3D> cage_points = SeedDetector::generate_cage(circle_vec, decimation);

  // Convert output: C++ vector of Points -> R DataFrame
  size_t total_points = cage_points.size();
  Rcpp::NumericVector out_X(total_points);
  Rcpp::NumericVector out_Y(total_points);
  Rcpp::NumericVector out_Z(total_points);

  for (size_t i = 0; i < total_points; ++i) {
    out_X[i] = cage_points[i].X;
    out_Y[i] = cage_points[i].Y;
    out_Z[i] = cage_points[i].Z;
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("X") = out_X,
    Rcpp::Named("Y") = out_Y,
    Rcpp::Named("Z") = out_Z
  );
}

Rcpp::DataFrame detect_tree_circles_cpp(Rcpp::DataFrame wood_df, double resolution = 0.05, int connectivity = 26, int num_ransac_iterations = 400, double inlier_threshold = 0.02, int min_cluster_size = 20)
{
  PointCloud wood(wood_df);

  std::vector<SeedDetectorGeometries::Circle> circles = SeedDetector::detect_tree_circles(
    wood,
    resolution,
    connectivity,
    num_ransac_iterations,
    inlier_threshold,
    min_cluster_size
  );

  size_t n_circles = circles.size();

  Rcpp::NumericVector X(n_circles);
  Rcpp::NumericVector Y(n_circles);
  Rcpp::NumericVector Z(n_circles);
  Rcpp::NumericVector R(n_circles);
  Rcpp::IntegerVector id(n_circles);

  for (size_t i = 0; i < n_circles; ++i)
  {
    X[i] = circles[i].X;
    Y[i] = circles[i].Y;
    Z[i] = circles[i].Z;
    R[i] = circles[i].R;
    id[i] = circles[i].id;
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("X") = X,
    Rcpp::Named("Y") = Y,
    Rcpp::Named("Z") = Z,
    Rcpp::Named("R") = R,
    Rcpp::Named("id") = id
  );
}


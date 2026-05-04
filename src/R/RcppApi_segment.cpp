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


// Type aliases for clarity
using GraphPtr = Rcpp::XPtr<arbor::segment::Graph>;
using GraphCache = std::pair<arbor::segment::Graph::DistanceVector, arbor::segment::Graph::PredecessorMap>;
using DF = Rcpp::DataFrame;

namespace arbor::segment
{
std::vector<int>  accumulate_passages(const PointCloud& core, const PointCloud& ground, const settings::GraphParameters& params);
std::vector<double>  dist2root(const PointCloud& core, const PointCloud& dtm, const settings::GraphParameters& params);
std::vector<bool> assign_wood_from_passage(const PointCloud& pc, const arbor::settings::SemanticParameters& params);
std::vector<bool> assign_wood_from_high_likelihood(const PointCloud& pc, const arbor::settings::SemanticParameters& params);
std::vector<bool> assign_wood_from_medium_likelihood(const PointCloud& pc, const arbor::settings::SemanticParameters& params);
std::vector<bool> assign_wood_from_wood_dilatation(const PointCloud& pc, const arbor::settings::SemanticParameters& params);
Graph* build_semantic_graph(const PointCloud& core, const PointCloud& target, const PointCloud& gnd, const arbor::settings::GraphParameters& params);
Graph* build_instance_graph(const PointCloud& core, const PointCloud& seeds, const arbor::settings::GraphParameters& params);
}

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
  for (size_t i = 0 ; i < p.size() ; i++) {
    if (p.get_treeid(i) == -1) {
      p.set_treeid(i, NA_INTEGER);
    }
  }
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

Rcpp::IntegerVector accumulate_passages_cpp(DF core, DF gnd, Rcpp::List params)
{
  arbor::settings::GraphParameters gparams = extract_pathfinder_params(params);
  PointCloud p(core);
  PointCloud s(gnd);
  std::vector<int> ans = arbor::segment::accumulate_passages(p, s, gparams);
  return Rcpp::IntegerVector(ans.begin(), ans.end());
}

Rcpp::NumericVector dist2root(DF core, DF gnd, Rcpp::List params)
{
  arbor::settings::GraphParameters gparams = extract_pathfinder_params(params);
  PointCloud p(core);
  PointCloud s(gnd);
  std::vector<double> ans = arbor::segment::dist2root(p, s, gparams);
  return Rcpp::NumericVector(ans.begin(), ans.end());
}


Rcpp::LogicalVector assign_wood_from_passage_cpp(DF core, Rcpp::List params)
{
  arbor::settings::SemanticParameters sparams = extract_semantic_params(params);
  PointCloud p(core);
  std::vector<bool> ans = arbor::segment::assign_wood_from_passage(p, sparams);
  return Rcpp::LogicalVector(ans.begin(), ans.end());
}

Rcpp::LogicalVector assign_wood_from_high_likelihood_cpp(DF core, Rcpp::List params)
{
  arbor::settings::SemanticParameters sparams = extract_semantic_params(params);
  PointCloud p(core);
  std::vector<bool> ans = arbor::segment::assign_wood_from_high_likelihood(p, sparams);
  return Rcpp::LogicalVector(ans.begin(), ans.end());
}

Rcpp::LogicalVector assign_wood_from_medium_likelihood_cpp(DF core, Rcpp::List params)
{
  arbor::settings::SemanticParameters sparams = extract_semantic_params(params);
  PointCloud p(core);
  std::vector<bool> ans = arbor::segment::assign_wood_from_medium_likelihood(p, sparams);
  return Rcpp::LogicalVector(ans.begin(), ans.end());
}

Rcpp::LogicalVector assign_wood_from_wood_dilatation_cpp(DF core, Rcpp::List params)
{
  arbor::settings::SemanticParameters sparams = extract_semantic_params(params);
  PointCloud p(core);
  std::vector<bool> ans = arbor::segment::assign_wood_from_wood_dilatation(p, sparams);
  return Rcpp::LogicalVector(ans.begin(), ans.end());
}

SEXP build_semantic_graph(DF dec, DF targets, DF gnd, Rcpp::List params)
{
  arbor::settings::GraphParameters p = extract_pathfinder_params(params);

  PointCloud core(dec);
  PointCloud trgt(targets);
  PointCloud ground(gnd);

  arbor::segment::Graph* graph = arbor::segment::build_semantic_graph(core, trgt, ground, p);
  GraphPtr ptr(graph, true);
  return ptr;
}

SEXP build_instance_graph(DF dec, DF seed, Rcpp::List params)
{
  arbor::settings::GraphParameters p = extract_pathfinder_params(params);

  PointCloud core(dec);
  PointCloud seeds(seed);

  arbor::segment::Graph* graph = arbor::segment::build_instance_graph(core, seeds, p);

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
  arbor::segment::Graph::GraphCache cache = graph->compute_distances(start_node);

  // Parallel loop over goal nodes
  #pragma omp parallel
  {
    std::vector<int> local_passage(num_points, 0);  // thread-local counts

    #pragma omp for schedule(dynamic, 100)
    for (int i = 0; i < n_goals; ++i)
    {
      arbor::segment::Graph::NodeId goal  = goal_nodes[i];

      auto [path, cost] = graph->findPath(start_node, goal, cache);

      for (size_t j = 0; j < path.size(); ++j)
      {
        arbor::segment::Graph::NodeId id = path[j];
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
  arbor::segment::Graph::NodeIDs node_ids(ids.begin(), ids.end());

  std::vector<double> distances;
  arbor::segment::Graph::NodeIDs closest_nodeids;

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

  std::vector<arbor::seeds::Circle> circle_vec;
  circle_vec.reserve(n);
  for (int i = 0; i < n; ++i) {

    circle_vec.push_back(arbor::seeds::Circle(X[i], Y[i], Z[i], R[i], id[i]));
  }

  // Call pure C++ function
  std::vector<arbor::seeds::Point3D> cage_points = arbor::seeds::SeedDetector::generate_cage(circle_vec, decimation);

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

  std::vector<arbor::seeds::Circle> circles = arbor::seeds::SeedDetector::detect_tree_circles(
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

#endif

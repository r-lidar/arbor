/*
 Pathfinding Logic Explanation:

 Dijkstra's Algorithm (compute_distances):
 The core of the pathfinding logic uses Dijkstra's algorithm, which computes the shortest path from
 a start node to all other nodes in the graph.
 It uses a priority queue (min-heap) to always expand the node with the smallest known distance.
 For each node, it checks all its neighbors and updates their distances if a shorter path is found.
 The algorithm continues until all nodes are processed or the queue is empty.

 Reconstructing the Path (findPath):
 After computing the shortest distances, the function findPath reconstructs the path from the start
 node to the goal node by following the predecessors.
 If a path exists, the nodes are added to the path in reverse order and then reversed at the end to
 obtain the correct order.

 Distance Matrix Calculation (getDistanceMatrix):
 This method computes the shortest paths for multiple start and goal nodes and returns a distance
 matrix, where each entry represents the shortest path between a start node and a goal node.
 It uses OpenMP to parallelize the calculation for efficiency when dealing with multiple nodes.
 */

#include <Rcpp.h>
#include "Graph.h"

// [[Rcpp::export]]
Rcpp::List findPaths(Rcpp::DataFrame graph_df, Rcpp::IntegerVector start_node_ids, Rcpp::IntegerVector goal_node_ids)
{
  Rcpp::IntegerVector from_nodes = graph_df["from"];
  Rcpp::IntegerVector to_nodes   = graph_df["to"];
  Rcpp::NumericVector edge_costs = graph_df["cost"];

  // Construct Graph using raw pointers (no data copy)
  Graph graph(from_nodes.begin(), to_nodes.begin(), edge_costs.begin(), from_nodes.size());

  // Precompute distances
  std::unordered_map<int, std::pair<std::vector<float>, std::unordered_map<int, int>>> precomputed_data;
  for (int start_node_id : start_node_ids)
  {
    if (precomputed_data.find(start_node_id) == precomputed_data.end())
    {
      precomputed_data[start_node_id] = graph.compute_distances(start_node_id);
    }
  }

  int num_pairs = start_node_ids.size();
  Rcpp::List paths(num_pairs);
  Rcpp::NumericVector total_costs(num_pairs);

  for (int i = 0; i < num_pairs; ++i)
  {
    auto result = graph.findPath(start_node_ids[i], goal_node_ids[i],
                                 precomputed_data[start_node_ids[i]]);
    paths[i] = result.first;
    total_costs[i] = result.second;
  }

  return Rcpp::List::create(
    Rcpp::_["paths"] = paths,
    Rcpp::_["total_costs"] = total_costs);
}

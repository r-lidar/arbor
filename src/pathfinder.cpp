#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <unordered_map>

#include "myomp.h" // Include OpenMP header for parallel processing

#include <Rcpp.h>

typedef std::vector<std::vector<float>> Matrix;

/*
 Pathfinding Logic Explanation:

 Dijkstra's Algorithm (compute_distances):
 The core of the pathfinding logic uses Dijkstra's algorithm, which computes the shortest path from a start node to all other nodes in the graph.
 It uses a priority queue (min-heap) to always expand the node with the smallest known distance.
 For each node, it checks all its neighbors and updates their distances if a shorter path is found.
 The algorithm continues until all nodes are processed or the queue is empty.

 Reconstructing the Path (findPath):
 After computing the shortest distances, the function findPath reconstructs the path from the start node to the goal node by following the predecessors.
 If a path exists, the nodes are added to the path in reverse order and then reversed at the end to obtain the correct order.

 Distance Matrix Calculation (getDistanceMatrix):
 This method computes the shortest paths for multiple start and goal nodes and returns a distance matrix, where each entry represents the shortest path between a start node and a goal node.
 It uses OpenMP to parallelize the calculation for efficiency when dealing with multiple nodes.
 */

struct Node {
  int destination_node_id;
  float edge_cost;
};


class Graph
{
public:
  std::vector<std::vector<Node>> adjacency_list;

  Graph(const Rcpp::DataFrame& graph_df)
  {
    Rcpp::IntegerVector from_nodes = graph_df["from"];
    Rcpp::IntegerVector to_nodes = graph_df["to"];
    Rcpp::NumericVector edge_costs = graph_df["cost"];

    int max_node = std::max(*std::max_element(from_nodes.begin(), from_nodes.end()), *std::max_element(to_nodes.begin(), to_nodes.end()));
    adjacency_list.resize(max_node + 1); // Resize the adjacency list to accommodate all nodes

    for (size_t i = 0; i < from_nodes.size(); ++i)
    {
      add_edge(from_nodes[i], to_nodes[i], edge_costs[i]);
    }
  }

  inline void add_edge(int source_node_id, int destination_node_id, float edge_cost)
  {
    adjacency_list[source_node_id].push_back({destination_node_id, edge_cost});
  }

  // Dijkstra's Algorithm to compute the shortest distances from the start node to all other nodes
  std::pair<std::vector<float>, std::unordered_map<int, int>> compute_distances(int start_node_id)
  {
    std::vector<float> distances(adjacency_list.size(), std::numeric_limits<float>::infinity()); // Initialize distances to infinity
    distances[start_node_id] = 0.0f;                                                             // Distance from the start node to itself is 0

    std::unordered_map<int, int> predecessors; // To store the previous node on the shortest path to reconstruct the path

    using DijkstraNode = std::pair<float, int>; // Pair of distance and node ID for the priority queue
    std::priority_queue<DijkstraNode, std::vector<DijkstraNode>, std::greater<>> open_set; // Min-heap priority queue
    open_set.push({0.0f, start_node_id}); // Push the start node with distance 0

    // Main loop of Dijkstra's Algorithm
    while (!open_set.empty())
    {
      auto [current_cost, current_node_id] = open_set.top(); // Get the node with the smallest distance
      open_set.pop();

      // If this cost is already greater than the stored distance, skip it
      if (current_cost > distances[current_node_id]) continue;

      // Check all neighbors of the current node
      for (const Node& neighbor : adjacency_list[current_node_id])
      {
        float new_cost = current_cost + neighbor.edge_cost; // Calculate the cost to reach this neighbor
        // If the new cost is smaller than the current recorded distance, update it
        if (new_cost < distances[neighbor.destination_node_id])
        {
          distances[neighbor.destination_node_id] = new_cost;
          predecessors[neighbor.destination_node_id] = current_node_id; // Record the predecessor to reconstruct the path
          open_set.push({new_cost, neighbor.destination_node_id}); // Add this neighbor to the open set
        }
      }
    }

    return {distances, predecessors}; // Return the distances and predecessors
  }

  // Method to find the shortest path from start_node_id to goal_node_id using precomputed distances
  std::pair<std::vector<int>, float> findPath(int start_node_id, int goal_node_id, const std::pair<std::vector<float>, std::unordered_map<int, int>>& precomputed_data)
  {
    const auto& [distances, predecessors] = precomputed_data; // Retrieve the precomputed distances and predecessors

    // If the goal node is unreachable, return an empty path and a cost of -1
    if (distances[goal_node_id] == std::numeric_limits<float>::infinity()) return {{}, -1.0f};

    // Reconstruct the path by following the predecessors from the goal node
    std::vector<int> path;
    for (int node = goal_node_id; node != start_node_id; node = predecessors.at(node))
    {
      path.push_back(node); // Add the node to the path
      if (predecessors.find(node) == predecessors.end()) return {{}, -1.0f}; // If no predecessor found, return empty (no path)
    }
    path.push_back(start_node_id);          // Add the start node to the path
    std::reverse(path.begin(), path.end()); // Reverse the path to get it from start to goal
    return {path, distances[goal_node_id]}; // Return the path and the total cost
  }

  // Method to get a distance matrix for multiple start and goal nodes
  std::vector<std::vector<float>> getDistanceMatrix(const std::vector<int>& start_node_ids, const std::vector<int>& goal_node_ids)
  {
    int num_start_nodes = start_node_ids.size();
    int num_goal_nodes = goal_node_ids.size();
    Matrix distance_matrix(num_start_nodes, std::vector<float>(num_goal_nodes, -1.0f)); // Initialize distance matrix with -1

    #pragma omp parallel for
    for (size_t i = 0; i < num_start_nodes; ++i)
    {
      auto [distances, _] = compute_distances(start_node_ids[i]); // Compute distances from the start node
      for (size_t j = 0; j < num_goal_nodes; ++j) {
        distance_matrix[i][j] = distances[goal_node_ids[j]]; // Store the distance to each goal node
      }
    }

    return distance_matrix; // Return the filled distance matrix
  }
};

// [[Rcpp::export]]
Rcpp::List findPaths(Rcpp::DataFrame graph_df, Rcpp::IntegerVector start_node_ids, Rcpp::IntegerVector goal_node_ids)
{
  Graph graph(graph_df); // Create the graph from the DataFrame

  // Precompute the distances and predecessors for all start nodes
  std::unordered_map<int, std::pair<std::vector<float>, std::unordered_map<int, int>>> precomputed_data;
  for (int start_node_id : start_node_ids)
  {
    if (precomputed_data.find(start_node_id) == precomputed_data.end())
    {
      precomputed_data[start_node_id] = graph.compute_distances(start_node_id); // Compute distances only once for each start node
    }
  }

  int num_pairs = start_node_ids.size();
  Rcpp::List paths(num_pairs);
  Rcpp::NumericVector total_costs(num_pairs);

  // For each pair of start and goal nodes, find the shortest path
  for (int i = 0; i < num_pairs; ++i)
  {
    auto result = graph.findPath(start_node_ids[i], goal_node_ids[i], precomputed_data[start_node_ids[i]]);
    paths[i] = result.first; // Store the path
    total_costs[i] = result.second; // Store the total cost of the path
  }

  return Rcpp::List::create(Rcpp::_["paths"] = paths, Rcpp::_["total_costs"] = total_costs); // Return the results
}

// [[Rcpp::export]]
Rcpp::NumericMatrix get_distance_matrix(Rcpp::DataFrame graph_df, Rcpp::IntegerVector start_node_ids, Rcpp::IntegerVector goal_node_ids)
{
  Graph graph(graph_df); // Create the graph from the DataFrame

  std::vector<int> start_nodes = Rcpp::as<std::vector<int>>(start_node_ids); // Convert Rcpp vector to std::vector
  std::vector<int> goal_nodes = Rcpp::as<std::vector<int>>(goal_node_ids); // Convert Rcpp vector to std::vector
  auto distance_matrix = graph.getDistanceMatrix(start_nodes, goal_nodes); // Compute the distance matrix

  Rcpp::NumericMatrix result(start_nodes.size(), goal_nodes.size()); // Prepare the result matrix

  // Fill the result matrix with the computed distances
  for (size_t i = 0; i < start_nodes.size(); ++i)
  {
    for (size_t j = 0; j < goal_nodes.size(); ++j)
    {
      result(i, j) = distance_matrix[i][j];
    }
  }

  // Set row and column names for the result matrix
  rownames(result) = Rcpp::as<Rcpp::CharacterVector>(start_node_ids);
  colnames(result) = Rcpp::as<Rcpp::CharacterVector>(goal_node_ids);

  return result; // Return the distance matrix
}

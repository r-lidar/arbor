#include "Graph.h"

Graph::Graph(const int* from, const int* to, const double* cost, size_t n_edges)
{
  if (n_edges == 0)
    return;

  // Compute the max node ID to size adjacency list
  int max_node = 0;
  for (size_t i = 0; i < n_edges; ++i)
  {
    if (from[i] > max_node) max_node = from[i];
    if (to[i] > max_node)   max_node = to[i];
  }

  adjacency_list.resize(max_node + 1);

  for (size_t i = 0; i < n_edges; ++i)
  {
    add_edge(from[i], to[i], static_cast<float>(cost[i]));
  }
}

inline void Graph::add_edge(int source_node_id, int destination_node_id, float edge_cost)
{
  adjacency_list[source_node_id].push_back({destination_node_id, edge_cost});
}

std::pair<std::vector<float>, std::unordered_map<int, int>>
  Graph::compute_distances(int start_node_id)
  {
    using namespace std::chrono;
    auto t0 = high_resolution_clock::now();

    std::vector<float> distances(adjacency_list.size(), std::numeric_limits<float>::infinity());
    distances[start_node_id] = 0.0f;

    std::unordered_map<int, int> predecessors;

    using DijkstraNode = std::pair<float, int>;
    std::priority_queue<DijkstraNode, std::vector<DijkstraNode>, std::greater<>> open_set;
    open_set.push({0.0f, start_node_id});

    while (!open_set.empty())
    {
      auto [current_cost, current_node_id] = open_set.top();
      open_set.pop();

      if (current_cost > distances[current_node_id]) continue;

      for (const Node& neighbor : adjacency_list[current_node_id])
      {
        float new_cost = current_cost + neighbor.edge_cost;
        if (new_cost < distances[neighbor.destination_node_id])
        {
          distances[neighbor.destination_node_id] = new_cost;
          predecessors[neighbor.destination_node_id] = current_node_id;
          open_set.push({new_cost, neighbor.destination_node_id});
        }
      }
    }

    auto t1 = high_resolution_clock::now();
    double elapsed_s = duration_cast<seconds>(t1 - t0).count();
    // Rcpp::Rcout << "[Dijkstra] Computation time: " << elapsed_s << " s" << std::endl;

    return {distances, predecessors};
  }

std::pair<std::vector<int>, float>
  Graph::findPath(int start_node_id, int goal_node_id,
                  const std::pair<std::vector<float>, std::unordered_map<int, int>>& precomputed_data)
  {
    const auto& [distances, predecessors] = precomputed_data;

    if (distances[goal_node_id] == std::numeric_limits<float>::infinity())
      return {{}, -1.0f};

    std::vector<int> path;
    int node = goal_node_id;

    while (node != start_node_id)
    {
      if (predecessors.find(node) == predecessors.end())
        return {{}, -1.0f};
      path.push_back(node);
      node = predecessors.at(node);
    }

    path.push_back(start_node_id);
    std::reverse(path.begin(), path.end());

    return {path, distances[goal_node_id]};
  }

std::vector<std::vector<float>>
  Graph::getDistanceMatrix(const std::vector<int>& start_node_ids, const std::vector<int>& goal_node_ids)
  {
    int num_start_nodes = start_node_ids.size();
    int num_goal_nodes = goal_node_ids.size();
    Matrix distance_matrix(num_start_nodes, std::vector<float>(num_goal_nodes, -1.0f));

    #pragma omp parallel for
    for (size_t i = 0; i < num_start_nodes; ++i)
    {
      auto [distances, _] = compute_distances(start_node_ids[i]);
      for (size_t j = 0; j < num_goal_nodes; ++j)
      {
        distance_matrix[i][j] = distances[goal_node_ids[j]];
      }
    }

    return distance_matrix;
  }

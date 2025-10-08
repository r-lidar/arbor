#include <limits>
#include <algorithm>

#include "Graph.h"

// --- Constructor ---
Graph::Graph(const int* from, const int* to, const double* cost, size_t n_edges)
{
  if (n_edges == 0) return;

  // Determine the maximum node ID to size adjacency list
  NodeId max_node = 0;
  for (size_t i = 0; i < n_edges; ++i)
    max_node = std::max({max_node, from[i], to[i]});

  adjacency_list.resize(max_node + 1);

  for (size_t i = 0; i < n_edges; ++i)
    add_edge(from[i], to[i], static_cast<EdgeCost>(cost[i]));
}

// --- Add edge ---
inline void Graph::add_edge(NodeId source, NodeId destination, EdgeCost cost)
{
  adjacency_list[source].push_back({destination, cost});
}

// --- Dijkstra’s algorithm ---
std::pair<DistanceVector, PredecessorMap>
  Graph::compute_distances(NodeId start) const
  {
    DistanceVector distances(adjacency_list.size(), std::numeric_limits<Cost>::infinity());
    distances[start] = 0.0f;

    PredecessorMap predecessors;

    using QueueNode = std::pair<Cost, NodeId>;
    std::priority_queue<QueueNode, std::vector<QueueNode>, std::greater<>> open_set;
    open_set.push({0.0f, start});

    while (!open_set.empty())
    {
      auto [current_cost, current] = open_set.top();
      open_set.pop();

      if (current_cost > distances[current]) continue;

      for (const Node& neighbor : adjacency_list[current])
      {
        Cost new_cost = current_cost + neighbor.cost;
        if (new_cost < distances[neighbor.destination])
        {
          distances[neighbor.destination] = new_cost;
          predecessors[neighbor.destination] = current;
          open_set.push({new_cost, neighbor.destination});
        }
      }
    }

    return {distances, predecessors};
  }

// --- Path reconstruction ---
std::pair<Path, Cost>
  Graph::findPath(NodeId start, NodeId goal,
                  const std::pair<DistanceVector, PredecessorMap>& precomputed_data) const
  {
    const auto& [distances, predecessors] = precomputed_data;

    if (distances[goal] == std::numeric_limits<Cost>::infinity())
      return {{}, -1.0f};

    Path path;
    for (NodeId node = goal; node != start; )
    {
      auto it = predecessors.find(node);
      if (it == predecessors.end())
        return {{}, -1.0f};
      path.push_back(node);
      node = it->second;
    }

    path.push_back(start);
    std::reverse(path.begin(), path.end());
    return {path, distances[goal]};
  }

// --- Compute matrix of distances ---
Matrix Graph::getDistanceMatrix(const std::vector<NodeId>& start_nodes,
                                const std::vector<NodeId>& goal_nodes) const
{
  const size_t num_start = start_nodes.size();
  const size_t num_goal = goal_nodes.size();
  Matrix distance_matrix(num_start, std::vector<Cost>(num_goal, -1.0f));

  #pragma omp parallel for
  for (size_t i = 0; i < num_start; ++i)
  {
    auto [distances, _] = compute_distances(start_nodes[i]);
    for (size_t j = 0; j < num_goal; ++j)
      distance_matrix[i][j] = distances[goal_nodes[j]];
  }

  return distance_matrix;
}


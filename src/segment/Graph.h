#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <unordered_map>
#include <queue>

#include "myomp.h"

using NodeId = int;
using Cost = float;
using EdgeCost = float;

struct Node
{
  NodeId destination;
  EdgeCost cost;
};

using AdjacencyList = std::vector<std::vector<Node>>;
using DistanceVector = std::vector<Cost>;
using PredecessorMap = std::unordered_map<NodeId, NodeId>;
using Matrix = std::vector<std::vector<Cost>>;
using Path = std::vector<NodeId>;

class Graph
{
public:
  AdjacencyList adjacency_list;

  Graph() = default;
  Graph(const int* from, const int* to, const double* cost, size_t n_edges);
  inline void add_edge(NodeId source, NodeId destination, EdgeCost cost);
  std::pair<DistanceVector, PredecessorMap> compute_distances(NodeId start) const;
  std::pair<Path, Cost> findPath(NodeId start, NodeId goal, const std::pair<DistanceVector, PredecessorMap>& precomputed_data) const;
  Matrix getDistanceMatrix(const std::vector<NodeId>& start_nodes, const std::vector<NodeId>& goal_nodes) const;

  // For each target point, get the closet start point
void shortest_paths_from_ground(
    const std::vector<NodeId>& ground_nodes,
    std::vector<double>& distances,
    std::vector<NodeId>& closest_ground
) const;
};

#endif // GRAPH_H

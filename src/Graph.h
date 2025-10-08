#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <unordered_map>
#include <queue>
#include <limits>
#include <algorithm>
#include <chrono>

#include "myomp.h"

typedef std::vector<std::vector<float>> Matrix;

struct Node
{
  int destination_node_id;
  float edge_cost;
};

class Graph
{
public:
  std::vector<std::vector<Node>> adjacency_list;

  Graph(const int* from, const int* to, const double* cost, size_t n_edges);

  inline void add_edge(int source_node_id, int destination_node_id, float edge_cost);

  std::pair<std::vector<float>, std::unordered_map<int, int>> compute_distances(int start_node_id);

  std::pair<std::vector<int>, float> findPath(
      int start_node_id,
      int goal_node_id,
      const std::pair<std::vector<float>, std::unordered_map<int, int>>& precomputed_data);

  std::vector<std::vector<float>> getDistanceMatrix(
      const std::vector<int>& start_node_ids,
      const std::vector<int>& goal_node_ids);
};

#endif

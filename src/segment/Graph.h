#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <unordered_map>

namespace arbor::segment {

class Graph
{
public:
  using NodeId = int;
  using Cost = float;
  using DistanceVector = std::vector<Cost>;
  using PredecessorMap = std::unordered_map<NodeId, NodeId>;
  using Path = std::vector<NodeId>;
  using NodeIDs = std::vector<NodeId>;
  using GraphCache = std::pair<DistanceVector, PredecessorMap>;

  struct Node
  {
    NodeId destination;
    Cost cost;
  };

  using AdjacencyList = std::vector<std::vector<Node>>;

public:
  AdjacencyList adjacency_list;

  Graph() = default;
  void resize(size_t n, int k);
  void add_edge(NodeId source, NodeId destination, Cost cost);
  std::pair<DistanceVector, PredecessorMap> compute_distances(NodeId start) const;
  std::pair<Path, Cost> findPath(NodeId start, NodeId goal, const std::pair<DistanceVector, PredecessorMap>& precomputed_data) const;
  void shortest_paths_from_node(const NodeIDs& ground_nodes, std::vector<double>& distances, NodeIDs& closest_ground) const;
};

}

#endif // GRAPH_H

#pragma once

#include <Rcpp.h>

#include "Adaptor.h"
#include "Graph.h"


// Type aliases for readability
using PointCloud = PointCloudAdaptor;
using NodeId = int;
using Cost = float;
using EdgeCost = float;
using ListOrNull = Rcpp::Nullable<Rcpp::List>;

using GraphPtr = Rcpp::XPtr<Graph>;

/**
 * @class GraphBuilder
 * @brief Helper class to progressively build a multi-layer graph from point clouds.
 *
 * Layers:
 *  - Point network (bidirectional KNN)
 *  - Target network (one-way: points → targets)
 *  - Ground network (one-way: ground → points)
 *  - Master seed network (one-way: master → ground)
 */
class GraphBuilder
{
public:
  int k = 20;
  float max_gap = 1.0f;
  float power = 3.0f;
  float wood2wood = 0.1;
  float leaf2leaf = 20;
  float wood2leaf = 100;
  bool downward = false;

private:
  Graph* graph;
  int offset_points = 0;
  int offset_targets = 0;
  int offset_ground = 0;
  int total_nodes = 0;
  std::vector<bool> wood;

public:

  GraphBuilder();
  ~GraphBuilder();
  Graph* get_graph();

  void add_core_layer(const PointCloud& dec);
  void add_target_layer(const PointCloud& dec, const PointCloud& target);
  void add_ground_layer(const PointCloud& dec,  const PointCloud& ground);
  void add_seed_layer(const PointCloud& dec,  const PointCloud& seeds);
  void add_master_seed_layer(const PointCloud& gnd, const PointCloud& master_seed);

  void set_wood(const std::vector<bool>& w);
};

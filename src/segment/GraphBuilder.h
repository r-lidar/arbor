#ifndef GRAPHBUILDER_H
#define GRAPHBUILDER_H

#include "Adaptor.h"
#include "Graph.h"
#include "api.h"

class GraphBuilder
{
public:
  GraphBuilder(const GraphBuilderParams& p);
  ~GraphBuilder();
  Graph* get_graph();

  void add_core_layer(const PointCloud& dec);
  void add_target_layer(const PointCloud& dec, const PointCloud& target);
  void add_seed_layer(const PointCloud& dec,  const PointCloud& seeds);
  void add_master_seed_layer();

  void set_wood(const std::vector<bool>& x);

  int get_num_cores() const;
  int get_num_targets() const;
  int get_num_seeds() const;
  int get_num_master() const;
  std::pair<int, int> get_range_core() const;
  std::pair<int, int> get_range_targets() const;
  std::pair<int, int> get_range_seed() const;
  std::pair<int, int> get_range_master() const;

private:
  Graph* graph;
  int offset_points = -1;
  int offset_targets = -1;
  int offset_seeds = -1;
  int offset_master = -1;
  int total_core_nodes = 0;
  int total_target_nodes = 0;
  int total_seed_nodes = 0;
  int total_master_nodes = 0;
  int total_nodes = 0;
  std::vector<bool> wood;
  bool graph_owner = true;
  GraphBuilderParams params;

  void set_angle_penalty(const std::vector<float>& x);
};

#endif

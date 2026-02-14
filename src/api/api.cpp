#include "api.h"
#include "GraphBuilder.h"

std::vector<bool> homogeneization(const PointCloud& pc, double res, bool hybrid = true);

NodeIDs segment_instance(const PointCloud& core, const PointCloud& seeds, const GraphBuilderParams& params)
{
  std::vector<bool> keep = homogeneization(core, 0.05, true);
  PointCloud dec = core.subset(keep);

  size_t num_raw_pts = core.size();
  size_t num_points  = dec.size();
  size_t num_trees   = seeds.size();

  Graph* graph = build_instance_graph(dec, seeds, params);

  NodeIDs seeds_ids(num_trees);
  std::iota(seeds_ids.begin(), seeds_ids.end(), num_points-1);

  std::vector<double> distances;
  NodeIDs closest_nodeids;
  graph->shortest_paths_from_node(seeds_ids, distances, closest_nodeids);

  return(closest_nodeids);
}

Graph* build_semantic_graph(const PointCloud& core, const PointCloud& targets, const PointCloud& ground, const GraphBuilderParams& params)
{
  GraphBuilder builder(params);

  builder.add_core_layer(core);
  builder.add_target_layer(core, targets);
  builder.add_seed_layer(core, ground);
  builder.add_master_seed_layer();

  return builder.get_graph();
}

Graph* build_instance_graph(const PointCloud& core, const PointCloud& seeds, const GraphBuilderParams& params)
{
  std::vector<bool> wood; wood.reserve(core.size());
  for (size_t i = 0; i < core.size(); ++i) wood.push_back(core.is_wood(i));

  GraphBuilder builder(params);
  builder.set_wood(wood);

  builder.add_core_layer(core);
  builder.add_seed_layer(core, seeds);
  builder.add_master_seed_layer();

  return builder.get_graph();
}

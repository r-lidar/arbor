#include "GraphBuilder.h"

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

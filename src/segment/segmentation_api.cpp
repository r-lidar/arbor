#include "GraphBuilder.h"

Graph* build_semantic_graph(PointCloud& core, PointCloud& targets, PointCloud& ground, PointCloud& master_seed, const GraphBuilderParams& params)
{
  GraphBuilder builder(params);

  builder.add_core_layer(core);
  builder.add_target_layer(core, targets);
  builder.add_seed_layer(core, ground);
  builder.add_master_seed_layer(ground, master_seed);

  return builder.get_graph();
}

Graph* build_instance_graph(PointCloud& core, PointCloud& seeds, PointCloud& master_seed, const GraphBuilderParams& params)
{
  std::vector<bool> wood; wood.reserve(core.size());
  for (size_t i = 0; i < core.size(); ++i) wood.push_back(core.get_foliage(i) == 0);

  GraphBuilder builder(params);
  builder.set_wood(wood);

  builder.add_core_layer(core);
  builder.add_seed_layer(core, seeds);
  builder.add_master_seed_layer(seeds, master_seed);

  return builder.get_graph();
}

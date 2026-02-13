#ifndef SEGMENTATION_API_H
#define SEGMENTATION_API_H

#include "Adaptor.h"
#include "GraphBuilder.h"

Graph* build_semantic_graph(const PointCloud& core, const PointCloud& target, const PointCloud& gnd, const PointCloud& master_seed, const GraphBuilderParams& params);
Graph* build_instance_graph(const PointCloud& core, const PointCloud& seeds, const PointCloud& master_seed, const GraphBuilderParams& params);

#endif

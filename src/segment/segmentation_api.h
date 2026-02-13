#ifndef SEGMENTATION_API_H
#define SEGMENTATION_API_H

#include "Adaptor.h"
#include "GraphBuilder.h"

Graph* build_semantic_graph(PointCloud& core, PointCloud& target, PointCloud& gnd, PointCloud& master_seed, const GraphBuilderParams& params);
Graph* build_instance_graph(PointCloud& core, PointCloud& seeds, PointCloud& master_seed, const GraphBuilderParams& params);

#endif

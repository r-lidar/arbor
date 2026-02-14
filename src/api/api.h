#ifndef API_H
#define API_H

#include "Adaptor.h"
#include "GraphBuilder.h"

Graph* build_semantic_graph(const PointCloud& core, const PointCloud& target, const PointCloud& gnd, const GraphBuilderParams& params);
Graph* build_instance_graph(const PointCloud& core, const PointCloud& seeds, const GraphBuilderParams& params);

#endif

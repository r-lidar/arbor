#ifndef API_H
#define API_H

#include "Adaptor.h"
#include "GraphBuilder.h"

using Logger = std::function<void(const std::string&)>;

std::vector<int> segment_instance(const PointCloud& core, const PointCloud& seeds, const GraphBuilderParams& params, const Logger& logger = [](const std::string&) {});
Graph* build_semantic_graph(const PointCloud& core, const PointCloud& target, const PointCloud& gnd, const GraphBuilderParams& params);
Graph* build_instance_graph(const PointCloud& core, const PointCloud& seeds, const GraphBuilderParams& params);

#endif

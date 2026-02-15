#ifndef API_H
#define API_H

#include "Adaptor.h"
#include "Parameters.h"
#include "GraphBuilder.h"

using Logger = std::function<void(const std::string&)>;

std::vector<int> segment_instance(const PointCloud& core, const PointCloud& seeds, const GraphBuilderParams& params, const Logger& logger = [](const std::string&) {});
std::vector<int> accumulate_passages(const PointCloud& core, const PointCloud& ground, const GraphBuilderParams& params, const Logger& logger = [](const std::string&) {});
std::vector<bool> assign_wood_from_passage(const PointCloud& pc, const SemanticParams& params, const Logger& logger = [](const std::string&) {});
std::vector<bool> assign_wood_from_high_likelihood(const PointCloud& pc, const SemanticParams& params, const Logger& logger = [](const std::string&) {});
std::vector<bool> assign_wood_from_medium_likelihood(const PointCloud& pc, const SemanticParams& params, const Logger& logger = [](const std::string&) {});
std::vector<bool> assign_wood_from_wood_dilatation(const PointCloud& pc, const SemanticParams& params, const Logger& logger = [](const std::string&) {});

Graph* build_semantic_graph(const PointCloud& core, const PointCloud& target, const PointCloud& gnd, const GraphBuilderParams& params);
Graph* build_instance_graph(const PointCloud& core, const PointCloud& seeds, const GraphBuilderParams& params);

#endif

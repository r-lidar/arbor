#include "api.h"
#include "myomp.h"
#include "nanoflann.h"
#include "GraphBuilder.h"

#include <chrono>
#include <sstream>
#include <iomanip>

using KDTree  = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 3>;
using index_t = nanoflann::KNNResultSet<double>::IndexType;

Graph* build_instance_graph(const PointCloud& core, const PointCloud& seeds, const GraphBuilderParams& params)
{
  if (core.size() == 0)     throw std::runtime_error("build_instance_graph: core point cloud is empty.");
  if (seeds.size() == 0)    throw std::runtime_error("build_instance_graph: seeds point cloud is empty.");
  if (!core.has_foliage())  throw std::runtime_error("build_instance_graph: core point cloud is missing required 'foliage' attribute.");
  if (!seeds.has_treeid())  throw std::runtime_error("build_instance_graph: seed point cloud is missing required 'treeid' attribute.");

  std::vector<bool> wood; wood.reserve(core.size());
  for (size_t i = 0; i < core.size(); ++i) wood.push_back(core.is_wood(i));

  GraphBuilder builder(params);
  builder.set_wood(wood);

  builder.add_core_layer(core);
  builder.add_seed_layer(core, seeds);
  builder.add_master_seed_layer();

  return builder.get_graph();
}

void segment_instance(PointCloud& core, const PointCloud& seeds, const GraphBuilderParams& params, const Logger& logger)
{
  if (core.size() == 0)     throw std::runtime_error("segment_instance: core point cloud is empty.");
  if (seeds.size() == 0)    throw std::runtime_error("segment_instance: seeds point cloud is empty.");
  if (!core.has_foliage())  throw std::runtime_error("segment_instance: core point cloud is missing required 'foliage' attribute.");
  if (!seeds.has_treeid())  throw std::runtime_error("segment_instance: seed point cloud is missing required 'treeid' attribute.");

  const auto t0 = std::chrono::steady_clock::now();

  logger("Instance segmentation start");
  logger("Decimating the point cloud... (1/4)");

  // Decimation
  std::vector<bool> keep = homogeneization(core, params.decimation, true);
  PointCloud dec = core.subset(keep);

  size_t num_raw_pts = core.size();
  size_t num_points  = dec.size();
  size_t num_seeds   = seeds.size();

  logger("Constructing the graph (2/4)");

  // Build graph
  Graph* graph = build_instance_graph(dec, seeds, params);

  if (graph == nullptr) throw std::runtime_error("segment_instance: Failed to build graph (null pointer returned).");

  // Indexes of the seeds in the graph
  Graph::NodeIDs seeds_ids(num_seeds);
  std::iota(seeds_ids.begin(), seeds_ids.end(), num_points);

  logger("Pathfinder (3/4)");

  // Retrieve the closest seed for each point
  std::vector<double> distances;
  Graph::NodeIDs closest_nodeids;
  graph->shortest_paths_from_node(seeds_ids, distances, closest_nodeids);

  if (closest_nodeids.size() != num_points+num_seeds+1) throw std::runtime_error("segment_instance: Pathfinding returned incomplete results.");

  // Remap seed index -> seed id for each decimated point
  Graph::NodeIDs treeID(num_points, -1);
  for (size_t i = 0 ; i < num_points ; i++)
  {
    Graph::NodeId id = closest_nodeids[i];
    if (id != -1) treeID[i] = static_cast<Graph::NodeId>(seeds.get_treeid(id-num_points));
  }

  logger("Assigning tree IDs to the dense point cloud (4/4)");

  // Expand seed id from dec point cloud to core point cloud
  KDTree tree(3, dec, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  tree.buildIndex();
  nanoflann::SearchParameters nanoparams;
  nanoparams.sorted = false;

  Graph::NodeIDs ans(num_raw_pts);

  #pragma omp parallel
  {
    std::vector<index_t> idx(1);
    std::vector<double> dist(1);
    double q[3];

    #pragma omp for schedule(static)
    for (size_t i = 0; i < num_raw_pts; ++i)
    {
      core.get_point(i, q);
      nanoflann::KNNResultSet<double> resultSet(1);
      resultSet.init(idx.data(), dist.data());
      tree.findNeighbors(resultSet, q, nanoparams);
      index_t id = idx[0];
      ans[i] = treeID[id];
    }
  }

  for (size_t i = 0 ; i < core.size() ; i++)  core.set_treeid(i, ans[i]);

  const auto t1 = std::chrono::steady_clock::now();
  const std::chrono::duration<double> elapsed = t1 - t0;

  std::ostringstream oss;
  oss << std::fixed << std::setprecision(1) << elapsed.count();

  logger("Instance segmentation completed in " + oss.str() + " s");

  delete graph;
}

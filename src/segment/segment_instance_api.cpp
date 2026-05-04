#include "arbor.h"
#include "myomp.h"
#include "nanoflann.h"
#include "Grid3D.h"
#include "GraphBuilder.h"

#include <chrono>
#include <numeric>
#include <sstream>
#include <iomanip>

using KDTree  = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 3>;
using index_t = nanoflann::KNNResultSet<double>::IndexType;

namespace arbor::segment {

Graph* build_instance_graph(const PointCloud& core, const PointCloud& seeds, const settings::GraphParameters& params)
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

void fix_small_isolated_low_clusters(PointCloud& las, double res = 0.05, int min_size = 200)
{
  const size_t n = las.size();

  // group indices by tree ID
  std::unordered_map<int, std::vector<int>> tree_to_indices;
  for (size_t i = 0; i < n; ++i)
  {
    if (!las.is_wood(i))        continue;
    if (las.get_hag(i) >= 3.0)  continue;
    int id = las.get_treeid(i);
    if (id < 0)                 continue; // NA from R or -1
    tree_to_indices[id].push_back(static_cast<int>(i));
  }

  // PER-TREE
  for (auto& [id, indices] : tree_to_indices)
  {
    // Build a sub-cloud with Z scaled down by 0.1 to flatten the search
    PointCloud sub(indices.size());
    for (size_t k = 0; k < indices.size(); ++k)
    {
      int i = indices[k];
      sub.set_x(k, las.get_x(i));
      sub.set_y(k, las.get_y(i));
      sub.set_z(k, las.get_z(i) * 0.1);
    }

    // Run connected components on the flattened sub-cloud
    Grid3D grid(sub, res);
    std::vector<int> cluster_ids = grid.connected_components(26);

    // Count points per cluster
    std::unordered_map<int, int> counts;
    for (int cid : cluster_ids) counts[cid]++;

    if (counts.size() <= 1) continue;

    // Find the largest cluster
    int best_cluster = -1, best_count = -1;
    for (auto& [cid, cnt] : counts)
    {
      if (cnt > best_count) { best_count = cnt; best_cluster = cid; }
    }

    // Reclassify non-largest cluster points as foliage
    for (size_t k = 0; k < indices.size(); ++k)
    {
      if (cluster_ids[k] != best_cluster)
        las.set_foliage(indices[k], 1);
    }
  }
}

void segment_instance(PointCloud& core, const PointCloud& seeds, const settings::ArborParameters& params)
{
  if (core.size() == 0)     throw std::runtime_error("segment_instance: point cloud is empty.");
  if (seeds.size() == 0)    throw std::runtime_error("segment_instance: seeds point cloud contains 0 seed.");
  if (!core.has_foliage())  throw std::runtime_error("segment_instance: point cloud is missing required 'foliage' attribute.");
  if (!seeds.has_treeid())  throw std::runtime_error("segment_instance: seed point cloud is missing required 'treeid' attribute.");

  ServiceLocator::logger()("Partitioning...");
  core.partition([&](size_t i) { return core.get_hag(i) > params.global.cut_above_ground; });

  const auto t0 = std::chrono::steady_clock::now();

  ServiceLocator::logger()("Instance segmentation start");
  ServiceLocator::logger()("Decimating the point cloud... (1/4)");

  // Decimation
  std::vector<bool> keep = arbor::utils::homogeneization(core, params.pathfinder.decimation, true);
  PointCloud dec = core.subset(keep);

  size_t num_raw_pts = core.size();
  size_t num_points  = dec.size();
  size_t num_seeds   = seeds.size();

  ServiceLocator::logger()("Constructing the graph (2/4)");

  // Build graph
  Graph* graph = build_instance_graph(dec, seeds, params.pathfinder);

  if (graph == nullptr) throw std::runtime_error("segment_instance: Failed to build graph (null pointer returned).");

  // Indexes of the seeds in the graph
  Graph::NodeIDs seeds_ids(num_seeds);
  std::iota(seeds_ids.begin(), seeds_ids.end(), num_points);

  ServiceLocator::logger()("Pathfinder (3/4)");

  // Retrieve the closest seed for each point
  std::vector<double> distances;
  Graph::NodeIDs closest_nodeids;
  graph->shortest_paths_from_node(seeds_ids, distances, closest_nodeids);
  delete graph;

  if (closest_nodeids.size() != num_points+num_seeds+1) throw std::runtime_error("segment_instance: Pathfinding returned incomplete results.");

  // Remap seed index -> seed id for each decimated point
  Graph::NodeIDs treeID(num_points, -1);
  for (size_t i = 0 ; i < num_points ; i++)
  {
    Graph::NodeId id = closest_nodeids[i];
    if (id != -1) treeID[i] = static_cast<Graph::NodeId>(seeds.get_treeid(id-num_points));
  }

  ServiceLocator::logger()("Assigning tree IDs to the dense point cloud (4/4)");

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

  ServiceLocator::logger()("Fix low isolated wood clusters");
  fix_small_isolated_low_clusters(core);

  ServiceLocator::logger()("Instance segmentation completed in " + oss.str() + " s");
}


std::vector<double> dist2root(const PointCloud& core, const PointCloud& dtm, const settings::GraphParameters& params)
{
  if (core.size() == 0) throw std::runtime_error("dist2root(): point cloud is empty.");
  if (dtm.size() == 0)  throw std::runtime_error("dist2root(): seeds point cloud is empty.");

  ServiceLocator::logger()("Constructing the graph");
  GraphBuilder builder(params);
  builder.add_core_layer(core);
  builder.add_seed_layer(core, dtm);
  builder.add_master_seed_layer();

  Graph* graph = builder.get_graph();
  if (graph == nullptr) throw std::runtime_error("dist2root(): failed to build graph (null pointer).");

  // Node ID layout: [0, num_points) = core, [num_points, num_points+num_gnd) = ground, master_id last
  const size_t num_points = core.size();
  const size_t num_gnd    = dtm.size();
  const Graph::NodeId master_id = static_cast<Graph::NodeId>(num_points + num_gnd);

  // Run Dijkstra from master seed
  ServiceLocator::logger()("Computing shortest paths to ground");
  const Graph::GraphCache cache = graph->compute_distances(master_id);
  const auto& [graph_distances, predecessors] = cache;

  // For each core point, walk the predecessor chain and accumulate Euclidean distance
  std::vector<double> euclidean_distance_to_root(num_points, -1.0);

  for (size_t i = 0; i < num_points; ++i)
  {
    if (graph_distances[i] == std::numeric_limits<Graph::Cost>::infinity())
      continue; // unreachable

    double total_euclidean = 0.0;
    Graph::NodeId current = static_cast<Graph::NodeId>(i);

    while (true)
    {
      auto it = predecessors.find(current);
      if (it == predecessors.end()) break; // reached master or disconnected

      Graph::NodeId prev = it->second;

      // Only accumulate Euclidean distance between two core nodes.
      // Seed and master nodes are virtual (no real 3D coords) and have 0-cost edges.
      if (prev < static_cast<Graph::NodeId>(num_points) && current < static_cast<Graph::NodeId>(num_points))
      {
        double c[3], p[3];
        core.get_point(current, c);
        core.get_point(prev, p);

        const double dx = c[0] - p[0];
        const double dy = c[1] - p[1];
        const double dz = c[2] - p[2];
        total_euclidean += std::sqrt(dx*dx + dy*dy + dz*dz);
      }

      current = prev;
    }

    euclidean_distance_to_root[i] = total_euclidean;
  }

  delete graph;
  return euclidean_distance_to_root;
}


}

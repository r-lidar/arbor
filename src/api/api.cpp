#include "api.h"
#include "nanoflann.h"
#include "GraphBuilder.h"
#include "Grid3D.h"

#include <chrono>
#include <sstream>
#include <iomanip>

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 3> KDTree;
using index_t = nanoflann::KNNResultSet<double>::IndexType;
using GraphCache = std::pair<DistanceVector, PredecessorMap>;

std::vector<bool> homogeneization(const PointCloud& pc, double res, bool hybrid = true);
std::vector<bool> sor(const PointCloud& pc, unsigned int k, double m, int ncpu);

std::vector<int> segment_instance(const PointCloud& core, const PointCloud& seeds, const GraphBuilderParams& params, const Logger& logger)
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
  NodeIDs seeds_ids(num_seeds);
  std::iota(seeds_ids.begin(), seeds_ids.end(), num_points);

  logger("Pathfinder (3/4)");

  // Retrieve the closest seed for each point
  std::vector<double> distances;
  NodeIDs closest_nodeids;
  graph->shortest_paths_from_node(seeds_ids, distances, closest_nodeids);

  if (closest_nodeids.size() != num_points+num_seeds+1) throw std::runtime_error("segment_instance: Pathfinding returned incomplete results.");

  // Remap seed index -> seed id for each decimated point
  NodeIDs treeID(num_points, -1);
  for (size_t i = 0 ; i < num_points ; i++)
  {
    NodeId id = closest_nodeids[i];
    if (id != -1) treeID[i] = static_cast<NodeId>(seeds.get_treeid(id-num_points));
  }

  logger("Assigning tree IDs to the dense point cloud (4/4)");

  // Expand seed id from dec point cloud to core point cloud
  KDTree tree(3, dec, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  tree.buildIndex();
  nanoflann::SearchParameters nanoparams;
  nanoparams.sorted = false;

  NodeIDs ans(num_raw_pts);

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

  const auto t1 = std::chrono::steady_clock::now();
  const std::chrono::duration<double> elapsed = t1 - t0;

  std::ostringstream oss;
  oss << std::fixed << std::setprecision(1) << elapsed.count();

  logger("Instance segmentation completed in " + oss.str() + " s");

  return(ans);
}

std::vector<int> accumulate_passages(const PointCloud& core, const PointCloud& ground, const GraphBuilderParams& params, const Logger& logger)
{
  if (core.size() == 0)   throw std::runtime_error("segment_instance: core point cloud is empty.");
  if (ground.size() == 0) throw std::runtime_error("segment_instance: seeds point cloud is empty.");

  logger("Decimating the point cloud (1/9)");

  // Decimation
  std::vector<bool> keep1 = homogeneization(core, params.decimation, true);
  PointCloud dec = core.subset(keep1, true);


  logger("Discretizing scene space (2/9)");

  std::vector<bool> keep2 = homogeneization(core, params.space_res, false);
  PointCloud targets = core.subset(keep2, true);

  logger("Constructing the graph (3/9)");

  // Build graph
  Graph* graph = build_semantic_graph(dec, targets, ground, params);

  if (graph == nullptr) throw std::runtime_error("segment_instance: Failed to build graph (null pointer returned).");

  logger("Accumulating passages (4/9)");

  size_t num_raw_points = core.size();
  size_t num_points = dec.size();
  size_t num_target = targets.size();
  size_t num_gnd    = ground.size();

  std::vector<int> target_ids(num_target);
  std::vector<int> ground_ids(num_gnd);
  int master_id = num_points + num_target + num_gnd;

  std::iota(target_ids.begin(), target_ids.end(), num_points);
  std::iota(ground_ids.begin(), ground_ids.end(), num_points + num_target);

  // Global count vector
  std::vector<int> passage(num_points, 0);

  // Precompute distances for fast access
  GraphCache cache = graph->compute_distances(master_id);

  // Parallel loop over goal nodes
  #pragma omp parallel
  {
    std::vector<int> local_passage(num_points, 0);  // thread-local counts

    #pragma omp for schedule(dynamic, 100)
    for (int i = 0; i < num_target; ++i)
    {
      NodeId goal  = target_ids[i];

      auto [path, cost] = graph->findPath(master_id, goal, cache);

      for (size_t j = 0; j < path.size(); ++j)
      {
        NodeId id = path[j];
        if (id >= 0 && id < num_points)
          local_passage[id] += 1;
      }
    }

    // Merge results into global passage safely
    #pragma omp critical
    {
      for (int i = 0; i < num_points; ++i)
        passage[i] += local_passage[i];
    }
  }

  // Transfer passage values from dec back to core
  // Points not in dec get value 0
  std::vector<int> core_passage(num_raw_points, 0);

  size_t dec_idx = 0;
  for (size_t i = 0; i < num_raw_points; ++i)
  {
    if (keep1[i])
    {
      core_passage[i] = passage[dec_idx];
      ++dec_idx;
    }
  }

  return core_passage;
}

std::vector<bool> assign_wood_from_passage(const PointCloud& pc, const SemanticParams& params, const Logger& logger)
{
  if (pc.size() == 0)     throw std::runtime_error("assign_wood_from_passage: point cloud is empty.");
  if (!pc.has_passage())  throw std::runtime_error("assign_wood_from_passage: point cloud is missing required 'passage' attribute.");

  logger("Pathfinder-based wood segmentation (5/9)");

  // Filter pseudo-skeleton: points with passage > min_passage
  std::vector<bool> skeleton_mask(pc.size(), false);
  for (size_t i = 0; i < pc.size(); ++i)
  {
    if (pc.get_passage(i) > params.min_passage)
      skeleton_mask[i] = true;
  }
  PointCloud passages = pc.subset(skeleton_mask, true);

  if (passages.size() == 0) throw std::runtime_error("assign_wood_from_passage: no passage points found (no points with passage > min_passage)");

  logger("  Building KDtree");

  // Spatial index of the skeleton
  KDTree tree(3, pc, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  tree.buildIndex();
  nanoflann::SearchParameters nanoparams;
  nanoparams.sorted = false;

  logger("  k-nn search");

  // Each point close enough from a passage is assigned wood
  std::vector<bool> is_wood(pc.size(), false);

  // Precompute squared distance threshold to avoid repeated multiplication or sqrt
  const double dist_threshold_sq = params.wood_assignation_dist * params.wood_assignation_dist;
  const int k = params.wood_assignation_k;

  #pragma omp parallel
  {
    std::vector<index_t> idx(k);
    std::vector<double> dist(k);
    double q[3];

    #pragma omp for schedule(static)
    for (size_t i = 0; i < passages.size(); ++i)
    {
      passages.get_point(i, q);
      nanoflann::KNNResultSet<double> resultSet(k);
      resultSet.init(idx.data(), dist.data());
      tree.findNeighbors(resultSet, q, nanoparams);

      for (int j = 0 ; j < 10 ; j++)
      {
        if (dist[j] < dist_threshold_sq)
          is_wood[idx[j]] = true;
      }
    }
  }

  return is_wood;
}

std::vector<bool> assign_wood_from_high_likelihood(const PointCloud& pc, const SemanticParams& params, const Logger& logger)
{
  if (pc.size() == 0)     throw std::runtime_error("assign_wood_from_high_likelihood: point cloud is empty.");
  if (!pc.has_foliage())  throw std::runtime_error("assign_wood_from_high_likelihood: point cloud is missing required 'foliage' attribute.");
  if (!pc.has_pwood())    throw std::runtime_error("assign_wood_from_high_likelihood: point cloud is missing required 'pwood' attribute.");

  logger("High likelihood based wood segmentation (6/9)");

  // Extract only high likelihood + already wood in previous step (assign_wood_from_passage)
  std::vector<bool> mask(pc.size(), false);
  for (size_t i = 0; i < pc.size(); ++i) {
    if (pc.is_wood(i) || pc.get_pwood(i) > params.high_pwood_threshold)
      mask[i] = true;
  }
  PointCloud wood = pc.subset(mask, true);

  logger("  Connected components computing");

  // Connected components to detect big clusters
  Grid3D grid(wood, params.connected_components_res);
  std::vector<int> cluster_ids = grid.connected_components(26);

  logger("  Connected components filtering");

  // Remove small clusters
  int max_id = *std::max_element(cluster_ids.begin(), cluster_ids.end());
  std::vector<std::size_t> counts(max_id + 1, 0);
  for (int id : cluster_ids) {
    if (id != 0)
      ++counts[id];
  }
  for (int& id : cluster_ids) {
    if (id != 0 && counts[id] < params.connected_components_min)
      id = 0;
  }

  // Assign original point cloud with wood/foliage
  std::vector<bool> is_wood(pc.size(), false);
  size_t j = 0;
  for (size_t i = 0; i < pc.size(); ++i)
  {
    if (mask[i]) {
      // In R, a point only survives if its cluster_id > 0
      // AND its count >= min.
      int cid = cluster_ids[j];
      if (cid > 0 && counts[cid] >= (size_t)params.connected_components_min) {
        is_wood[i] = true;
      }

      ++j;
    }
  }

  return is_wood;
}

std::vector<bool> assign_wood_from_medium_likelihood(const PointCloud& pc, const SemanticParams& params, const Logger& logger)
{
  if (pc.size() == 0)     throw std::runtime_error("assign_wood_from_medium_likelihood: point cloud is empty.");
  if (!pc.has_foliage())  throw std::runtime_error("assign_wood_from_medium_likelihood: point cloud is missing required 'foliage' attribute.");
  if (!pc.has_pwood())    throw std::runtime_error("assign_wood_from_medium_likelihood: point cloud is missing required 'pwood' attribute.");

  logger("Medium likelihood based wood segmentation (7/9)");

  // Extract only medium likelihood + already wood in previous steps
  std::vector<bool> mask(pc.size(), false);
  for (size_t i = 0; i < pc.size(); ++i) {
    if (pc.is_wood(i) || (pc.get_pwood(i) > params.medium_pwood_threshold && pc.get_pwood(i) < params.high_pwood_threshold))
      mask[i] = true;
  }
  PointCloud wood = pc.subset(mask, true);

  logger("  sor noise segmentation");

  // SOR. Detect noise
  std::vector<bool> is_noise = sor(wood, params.medium_pwood_sor_k, params.medium_pwood_sor_m, 12);

  // Remove noise
  is_noise.flip();
  wood = wood.subset(is_noise, true);
  is_noise.flip();

  logger("  Connected components computing");

  // Connected components to detect big clusters
  Grid3D grid(wood, params.connected_components_res);
  std::vector<int> cluster_ids = grid.connected_components(26);

  logger("  Connected components filtering");

  // Remove small clusters
  int max_id = *std::max_element(cluster_ids.begin(), cluster_ids.end());
  std::vector<std::size_t> counts(max_id + 1, 0);
  for (int id : cluster_ids) {
    if (id != 0)
      ++counts[id];
  }
  for (int& id : cluster_ids) {
    if (id != 0 && counts[id] < params.connected_components_min)
      id = 0;
  }

  // Assign original point cloud with wood/foliage
  std::vector<bool> is_wood(pc.size(), false);
  size_t j = 0; // Index for wood_subset (mask)
  size_t k = 0; // Index for clustered_subset (survival_mask)

  for (size_t i = 0; i < pc.size(); ++i)
  {
    if (mask[i])
    {
      // Point was in the first subset. Did it survive SOR?
      if (!is_noise[j])
      {
        // Point was in the second subset. Did it survive Cluster Filtering?
        int cid = cluster_ids[k];
        if (cid > 0 && counts[cid] >= (size_t)params.connected_components_min) {
          is_wood[i] = true;
        }
        k++; // Increment k only if the point survived SOR
      }
      j++; // Increment j every time mask[i] is true
    }
  }

  return is_wood;
}

std::vector<bool> assign_wood_from_wood_dilatation(const PointCloud& pc, const SemanticParams& params, const Logger& logger)
{
  if (pc.size() == 0)     throw std::runtime_error("assign_wood_from_medium_likelihood: point cloud is empty.");
  if (!pc.has_foliage())  throw std::runtime_error("assign_wood_from_medium_likelihood: point cloud is missing required 'foliage' attribute.");

  // We look at the neighboring points of the wood.  Points close to the wood
  // are wood points too. This assigns extra wood point is the branches and remove
  // some false negatives

  logger("Dilatation based wood segmentation... (8/9)");

  // Extract wood points
  std::vector<bool> is_wood(pc.size(), false);
  for (size_t i = 0; i < pc.size(); ++i) {
    if (pc.is_wood(i))
      is_wood[i] = true;
  }
  PointCloud wood = pc.subset(is_wood, true);


  logger("  Building KDtree");

  // Build KDTree on wood points
  KDTree tree(3, pc, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  tree.buildIndex();
  nanoflann::SearchParameters nanoparams;
  nanoparams.sorted = false;


  logger("  knn search");

  const double max_dist_sq = params.wood_extra_reasignation_dist * params.wood_extra_reasignation_dist;
  const int k = params.wood_extra_reasignation_k;

  #pragma omp parallel
  {
    std::vector<index_t> idx(k);
    std::vector<double> dist(k);
    double q[3];

    #pragma omp for schedule(static)
    for (size_t i = 0; i < wood.size(); ++i)
    {
      wood.get_point(i, q);
      nanoflann::KNNResultSet<double> resultSet(k);
      resultSet.init(idx.data(), dist.data());
      tree.findNeighbors(resultSet, q, nanoparams);

      // If nearest wood point is within threshold distance, assign as wood
      for (int j = 0 ; j < 10 ; j++)
      if (dist[j] <= max_dist_sq)
        is_wood[idx[j]] = true;
    }
  }

  return is_wood;
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

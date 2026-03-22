#include "GraphBuilder.h"
#include "nanoflann.h"
#include "PointCloud.h"
#include "myomp.h"

#include <vector>
#include <unordered_map>
#include <cmath>

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>,PointCloud, 3>;

namespace arbor::segment {

// --------------------------------------------------------
// Constructor / Destructor
// --------------------------------------------------------

GraphBuilder::GraphBuilder(const settings::GraphParameters& p)
{
  set_angle_penalty(p.angle_penalty);
  params = p;
  graph = new Graph();
}

GraphBuilder::~GraphBuilder() { if (graph_owner) delete graph; }
Graph* GraphBuilder::get_graph() { graph_owner = false; return graph; }
void GraphBuilder::set_wood(const std::vector<bool>& x) { wood = x; }
void GraphBuilder::set_angle_penalty(const std::vector<float>& x)
{
  constexpr std::size_t expected_size = 181;

  if (x.size() != expected_size)
  {
    throw std::runtime_error(
        "Invalid angle penalty factor vector size: expected " +
          std::to_string(expected_size) +
          ", got " +
          std::to_string(x.size()) + "."
    );
  }

  params.angle_penalty = x;
}


// ---------------------------------------------------------
// 1. Core Layer (bidirectional)
// ---------------------------------------------------------

void GraphBuilder::add_core_layer(const PointCloud& core)
{
  if (total_core_nodes > 0)   throw std::runtime_error("Core layer already populated");
  if (total_target_nodes > 0) throw std::runtime_error("Core layer must be populated first");
  if (total_seed_nodes > 0)   throw std::runtime_error("Core layer must be populated first");
  if (total_master_nodes > 0) throw std::runtime_error("Core layer must be populated first");


  // Because self point is included in knn
  params.k++;

  int n_points = core.size();
  bool use_wood = wood.size() > 0;

  offset_points = 0;
  total_core_nodes = n_points;
  total_nodes = n_points;

  graph->ensure_size(total_nodes);

  // Build the KD-tree index
  KDTree index(3, core, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  // Num. threads used
  int n_threads = omp_get_max_threads();

  #pragma omp parallel
  {
    std::vector<size_t> idx(params.k);
    std::vector<double> dist(params.k);

    // Thread-local storage for edges
    std::vector<std::tuple<Graph::NodeId, Graph::NodeId, Graph::Cost>> local_edges;

    // For each point we connect to its knn. The current point is 'from'
    #pragma omp for schedule(static)
    for (int from = 0; from < n_points; ++from)
    {
      // Get the knn
      nanoflann::KNNResultSet<double> result(params.k);
      result.init(&idx[0], &dist[0]);
      double q[3];
      core.get_point(from, q);
      index.findNeighbors(result, q, nanoflann::SearchParameters());

      // For each knn, compute the cost to connect 'from' and 'to'
      for (int j = 0; j < params.k; ++j)
      {
        int to = idx[j];                        // This is our target index
        if (to == from) continue;               // If 'from' == 'to', skip because this is the 0-nn
        float cost = std::sqrt(dist[j]);        // The cost is the euclidean distance
        if (cost > params.max_gap) continue;    // If the cost is above a threshold; no connection
        cost = std::pow(cost, params.power);    // The cost is the cube of the eucliandian distance

        double coord_from[3];
        double coord_to[3];
        core.get_point(from, coord_from);
        core.get_point(to, coord_to);

        // Apply a extra cost factor base on the direction of the link
        // Moving upward is cheap. Downward is expensive.
        float dx = coord_from[0] - coord_to[0];
        float dy = coord_from[1] - coord_to[1];
        float dz = coord_from[2] - coord_to[2];
        float magnitude = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (magnitude < 1e-12) continue;
        float cos_theta = -dz / magnitude;
        if (params.downward) cos_theta = -cos_theta;
        float angle_deg = std::acos(std::clamp(cos_theta, -1.0f, 1.0f)) * 180.0f / M_PI;
        int angle = std::round(angle_deg);
        cost *= params.angle_penalty[angle_deg];

        // If we have a wood/foliage classification we apply extra cost factors
        if (use_wood)
        {
          bool is_wood1 = wood[from];
          bool is_wood2 = wood[to];
          if (is_wood1 && is_wood2) cost *= params.wood2wood;
          else if (!is_wood1 && !is_wood2) cost *= params.leaf2leaf;
          else if (is_wood1 && !is_wood2) cost *= params.wood2leaf;
        }

        // Add an edge per thread
        local_edges.emplace_back(from, to, cost);
      }
    }

    // Parallel reduction
    #pragma omp critical
    {
      for (auto& e : local_edges)
        graph->add_edge(std::get<0>(e), std::get<1>(e), std::get<2>(e));
    }
  }
}

// ---------------------------------------------------------
// 2. Target Layer (point → target)
// ---------------------------------------------------------

// Each target point is connected to its 1-nn core point

void GraphBuilder::add_target_layer(const PointCloud& core, const PointCloud& target)
{
  if (total_target_nodes > 0) throw std::runtime_error("Target layer already populated");
  if (total_core_nodes == 0)  throw std::runtime_error("Target layer must be populated after core layer");
  if (total_seed_nodes > 0)   throw std::runtime_error("Target layer must be populated before seed layer");
  if (total_master_nodes > 0) throw std::runtime_error("Target layer must be populated before master layer");

  int n_points = core.size();
  int n_target = target.size();

  offset_targets = total_nodes;
  total_target_nodes = n_points;
  total_nodes += n_target;

  graph->ensure_size(total_nodes);

  KDTree index(3, core, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  std::vector<size_t> idx(params.k);
  std::vector<double> dist(params.k);

  // For each target point, search the closest core point. The connection is
  // core point to target point
  for (int i = 0; i < n_target; ++i)
  {
    double q[3];
    target.get_point(i, q);
    nanoflann::KNNResultSet<double> result(params.k);
    result.init(&idx[0], &dist[0]);
    index.findNeighbors(result, q, nanoflann::SearchParameters());

    // Connection between point 'from' in the core network
    // to 'to' which is i + an offset to account for the existing
    // network. The cost is 0;
    int from = idx[0];
    int to = i + offset_targets;
    graph->add_edge(from, to, 0.0f);
  }
}

// ---------------------------------------------------------
// 3. Ground or seed Layer (ground → point)
// ---------------------------------------------------------

void GraphBuilder::add_seed_layer(const PointCloud& core, const PointCloud& seeds)
{
  if (total_seed_nodes > 0)    throw std::runtime_error("Seed layer already populated");
  if (total_core_nodes == 0)   throw std::runtime_error("Seed layer must be populated after core layer");
  if (total_master_nodes > 0)  throw std::runtime_error("Seed layer must be populated before master layer");

  int k = params.k_seed;
  int n_points = seeds.size();

  offset_seeds = total_nodes;
  total_seed_nodes = n_points;
  total_nodes += n_points;

  graph->ensure_size(total_nodes);

  KDTree index(3, core, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  std::vector<size_t> idx(k);
  std::vector<double> dist(k);

  for (int i = 0; i < n_points; ++i)
  {
    double q[3];
    seeds.get_point(i, q);
    nanoflann::KNNResultSet<double> result(k);
    result.init(&idx[0], &dist[0]);
    index.findNeighbors(result, q, nanoflann::SearchParameters());

    for (int j = 0; j < k; ++j)
    {
      float cost = std::sqrt(dist[j]);
      cost = std::pow(cost, params.power);    // The cost is the cube of the eucliandian distance
      int from = i + offset_seeds;
      int to = idx[j];
      graph->add_edge(from, to, cost);
    }
  }
}

// ---------------------------------------------------------
// 4. Master Seed Layer (master → all ground)
// ---------------------------------------------------------

// A master seed is connected to all ground points with cost 0

void GraphBuilder::add_master_seed_layer()
{
  if (total_master_nodes > 0)  throw std::runtime_error("Master layer already populated");
  if (total_core_nodes == 0)   throw std::runtime_error("Master layer must be populated after core layer");
  if (total_seed_nodes == 0)   throw std::runtime_error("Seed layer must be populated before master layer");

  offset_master = total_nodes;
  total_master_nodes = 1;
  total_nodes += 1;

  graph->ensure_size(total_nodes);

  for (int i = 0; i < total_seed_nodes; ++i)
  {
    graph->add_edge(offset_master, offset_seeds + i, 0.0f);
  }
}


int GraphBuilder::get_num_cores() const { return total_core_nodes; }
int GraphBuilder::get_num_targets() const { return total_target_nodes; }
int GraphBuilder::get_num_seeds() const { return total_seed_nodes; }
int GraphBuilder::get_num_master() const { return total_master_nodes; }
std::pair<int, int> GraphBuilder::get_range_core() const { return {offset_points, offset_points + total_core_nodes - 1}; }
std::pair<int, int> GraphBuilder::get_range_targets() const { return {offset_targets, offset_targets + total_target_nodes - 1}; }
std::pair<int, int> GraphBuilder::get_range_seed() const { return {offset_seeds, offset_seeds + total_seed_nodes - 1}; }
std::pair<int, int> GraphBuilder::get_range_master() const { return {offset_master, offset_master + total_master_nodes - 1 }; }

}

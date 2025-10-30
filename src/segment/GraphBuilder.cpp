#include "GraphBuilder.h"
#include "nanoflann.h"
#include "myomp.h"

#include <vector>
#include <unordered_map>
#include <cmath>

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>,PointCloud, 3>;

// --------------------------------------------------------
// Constructor / Destructor
// --------------------------------------------------------

GraphBuilder::GraphBuilder() { graph = new Graph(); }
GraphBuilder::~GraphBuilder() { if (graph_owner) delete graph; }
Graph* GraphBuilder::get_graph() { graph_owner = false; return graph; }
void GraphBuilder::set_wood(const std::vector<bool>& w) { wood = w; }

// ---------------------------------------------------------
// 1. Core Layer (bidirectional)
// ---------------------------------------------------------

void GraphBuilder::add_core_layer(const PointCloud& dec)
{
  // Because self point is included in knn
  k++;

  int n_points = dec.point_count();
  bool use_wood = wood.size() > 0;

  offset_points = 0;
  total_nodes = n_points;

  graph->ensure_size(total_nodes);

  // Angle cost factor lambda
  auto angle_penalty_factor = [](double angle_deg)
  {
    // precomputed std::log(100)/100 such a angle_penalty_factor(100) = 100
    constexpr double f = 0.046051;
    double y = std::exp(f * angle_deg);
    if (angle_deg > 100.0) y = 100.0;
    return y;
  };

  // Build the KD-tree index
  KDTree index(3, dec, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  // Num. threads used
  int n_threads = omp_get_max_threads();

  #pragma omp parallel
  {
    std::vector<size_t> idx(k);
    std::vector<double> dist(k);

    // Thread-local storage for edges
    std::vector<std::tuple<NodeId, NodeId, Cost>> local_edges;

    // For each point we connect to its knn. The current point is 'from'
    #pragma omp for schedule(static)
    for (int from = 0; from < n_points; ++from)
    {
      // Get the knn
      nanoflann::KNNResultSet<double> result(k);
      result.init(&idx[0], &dist[0]);
      double q[3];
      dec.get_point(from, q);
      index.findNeighbors(result, q, nanoflann::SearchParameters());

      // For each knn, compute the cost to connect 'from' and 'to'
      for (int j = 0; j < k; ++j)
      {
        int to = idx[j];                 // This is our target index
        if (to == from) continue;        // If 'from' == 'to', skip because this is the 0-nn
        float cost = std::sqrt(dist[j]); // The cost is the euclidean distance
        if (cost > max_gap) continue;    // If the cost is above a threshold; no connection
        cost = std::pow(cost, power);    // The cost is the cube of the eucliandian distance

        double coord_from[3];
        double coord_to[3];
        dec.get_point(from, coord_from);
        dec.get_point(to, coord_to);

        // Apply a extra cost factor base on the direction of the link
        // Moving upward is cheap. Downward is expensive.
        float dx = coord_from[0] - coord_to[0];
        float dy = coord_from[1] - coord_to[1];
        float dz = coord_from[2] - coord_to[2];
        float magnitude = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (magnitude < 1e-12) continue;
        float cos_theta = -dz / magnitude;
        if (downward) cos_theta = -cos_theta;
        float angle_deg = std::acos(std::clamp(cos_theta, -1.0f, 1.0f)) * 180.0f / M_PI;
        cost *= angle_penalty_factor(angle_deg);

        // If we have a wood/foliage classification we apply extra cost factors
        if (use_wood)
        {
          bool is_wood1 = wood[from];
          bool is_wood2 = wood[to];
          if (is_wood1 && is_wood2) cost *= wood2wood;
          else if (!is_wood1 && !is_wood2) cost *= leaf2leaf;
          else if (is_wood1 && !is_wood2) cost *= wood2leaf;
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

void GraphBuilder::add_target_layer(const PointCloud& dec, const PointCloud& target)
{
  int n_points = dec.point_count();
  int n_target = target.point_count();

  offset_targets = total_nodes;
  total_nodes += n_target;

  graph->ensure_size(total_nodes);

  KDTree index(3, dec, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  std::vector<size_t> idx(k);
  std::vector<double> dist(k);

  // For each target point, search the closest core point. The connection is
  // core point to target point
  for (int i = 0; i < n_target; ++i)
  {
    double q[3];
    target.get_point(i, q);
    nanoflann::KNNResultSet<double> result(k);
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

// Each ground point is connected to its k-nn core points

void GraphBuilder::add_ground_layer(const PointCloud& dec, const PointCloud& ground)
{
  int k = this->k*10;
  int n_gnd = ground.point_count();

  offset_ground = total_nodes;
  total_nodes += n_gnd;

  graph->ensure_size(total_nodes);

  KDTree index(3, dec, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  std::vector<size_t> idx(k);
  std::vector<double> dist(k);

  for (int i = 0; i < n_gnd; ++i)
  {
    double q[3];
    ground.get_point(i, q);
    nanoflann::KNNResultSet<double> result(k);
    result.init(&idx[0], &dist[0]);
    index.findNeighbors(result, q, nanoflann::SearchParameters());

    for (int j = 0; j < k; ++j)
    {
      float cost = std::sqrt(dist[j]);
      int from = i + offset_ground;
      int to = idx[j];
      graph->add_edge(from, to, cost);
    }
  }
}


void GraphBuilder::add_seed_layer(const PointCloud& dec, const PointCloud& seeds)
{
  int k = this->k*10;
  int n_gnd = seeds.point_count();

  offset_ground = total_nodes;
  total_nodes += n_gnd;

  graph->ensure_size(total_nodes);

  KDTree index(3, dec, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  std::vector<size_t> idx(k);
  std::vector<double> dist(k);

  for (int i = 0; i < n_gnd; ++i)
  {
    double q[3];
    seeds.get_point(i, q);
    nanoflann::KNNResultSet<double> result(k);
    result.init(&idx[0], &dist[0]);
    index.findNeighbors(result, q, nanoflann::SearchParameters());

    for (int j = 0; j < k; ++j)
    {
      float cost = std::sqrt(dist[j]);
      int from = i + offset_ground;
      int to = idx[j];
      graph->add_edge(from, to, cost);
    }
  }
}

// ---------------------------------------------------------
// 4. Master Seed Layer (master → all ground)
// ---------------------------------------------------------

// A master seed is connected to all ground points with cost 0

void GraphBuilder::add_master_seed_layer(const PointCloud& gnd, const PointCloud& master_seed)
{
  int n_gnd = gnd.point_count();

  int offset_seed = total_nodes;
  total_nodes += 1;

  graph->ensure_size(total_nodes);

  for (int i = 0; i < n_gnd; ++i)
  {
    graph->add_edge(offset_seed, offset_ground + i, 0.0f);
  }
}


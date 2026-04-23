#include "QSMbuilder.h"

#include <limits>
#include <vector>
#include <unordered_set>

namespace arbor::qsm {

void QSMbuilder::build(const PointCloud& tree)
{
  std::size_t n = tree.size();

  // Find the vertical bounds of the tree
  double min_z = std::numeric_limits<double>::max();
  double max_z = -std::numeric_limits<double>::max();

  for (std::size_t i = 0; i < n; ++i)
  {
    double current_z = tree.get_z(i);
    if (current_z < min_z) min_z = current_z;
    if (current_z > max_z) max_z = current_z;
  }

  // Calculate height and 1% threshold
  double tree_height = max_z - min_z;
  double trim_offset = tree_height * 0.01;
  double z_threshold = min_z + trim_offset;

  // Filter wood only and apply the dynamic trim
  std::vector<bool> wood_mask(n, false);
  for (std::size_t i = 0; i < n; ++i)
  {
    wood_mask[i] = tree.is_wood(i) && (tree.get_z(i) >= z_threshold);
  }

  PointCloud wood = tree.subset(wood_mask);
  n = wood.size();

  // Determine the geographic coordinates minimum
  // and center on 0,0,0 for numerical stability
  size_t min_idx = 0;
  min_z = tree.get_z(0);
  for (size_t i = 0; i < n; ++i)
  {
    double current_z = wood.get_z(i);
    if (current_z < min_z)
    {
      min_z = current_z;
      min_idx = i;
    }
  }
  double tx = wood.get_x(min_idx);
  double ty = wood.get_y(min_idx);
  double tz = wood.get_z(min_idx);
  wood.translate(tx, ty, tz);

  // Sometime the very bottom have two clusters of wood
  // this creates trouble. One of the two is removed to ensure
  // a single entry point for the QSM
  wood = clean_tree_butt(wood);
  n = wood.size();

  // REMOVED: it was not a good idea. This creates very bad trees ending and break in trunks
  // 3D smoothing for a better architecture computation
  // (but measure on the original on)
  //ServiceLocator::logger()("Smoothing 3D for better skeleton");
  //PointCloud swood = utils::smooth3d(wood, 0.04, 1);
  PointCloud& swood = wood;

  // Inspired by aRchi and needed to build the skeleton
  auto layers = this->layers(swood, params.qsm.step);
  auto clusters = this->clusters(swood, layers, params.qsm.cl_dist);

  if (layers.size() != clusters.size()) throw std::runtime_error("Internal error in QSMbuilder::build. layers and clusters have different sizes. Please report to info@r-lidar.com");
  if (layers.size() != wood.size()) throw std::runtime_error("Internal error in QSMbuilder::build. wood and layer have different sizes. Please report to info@r-lidar.com");

  // Convert std::vector<pair>
  std::vector<std::pair<int, int>> iter_cluster;
  iter_cluster.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    iter_cluster.emplace_back(layers[i].first, clusters[i].first);
  }

  // Build the QSM nodes
  build_skeleton(swood, iter_cluster, params.qsm.max_d);

  // Extremely rare case with so few points that we have no cluster
  // (seen once with a very bad DTM in Murray's data)
  if (graph.edges().size() == 0)
  {
    shift(tx, ty, tz);
    return;
  }

  // Connect the QSM nodes (sets parent_ID in each edge)
  compute_topology();

  // Fix root issue (rare)
  int n_root = count_nodes_connected_to_root();
  if (n_root == 0) throw std::runtime_error("Internal error in QSMbuilder::build. 0 root for this QSM. Please report to info@r-lidar.com.");
  if (n_root > 1)
  {
    ServiceLocator::logger()("Multiple nodes connected to root detected");
    fix_multiple_root();
  }

  compute_architecture(false);
  smooth_skeleton(params.qsm.smooth_steps, params.qsm.smooth_lambda, params.qsm.smooth_mu);
  detect_weird_butt();
  construct_radii(wood, params.qsm.apex);
  refine_radii(wood);
  estimate_prolongation(wood);
  prolongate(prolongation_distance);
  smooth_radii(15, 0.5, -0.7);
  shift(tx, ty, tz);
}

void QSMbuilder::shift(double tx, double ty, double tz)
{
  // Shift all node positions (a single update per node covers all incident edges)
  for (auto& [nid, ndata] : graph.nodes())
  {
    ndata.x += tx;
    ndata.y += ty;
    ndata.z += tz;
  }
}

int QSMbuilder::count_nodes_connected_to_root() const
{
  // Iterate through the public nodes map
  for (const auto& [id, data] : graph.nodes())
  {
    // The root is defined as a node with zero incoming edges
    if (graph.incoming_edges(id).empty())
    {
      // Return the count of its outgoing edges
      return graph.outgoing_edges(id).size();
    }
  }

  return 0;
}

int QSMbuilder::find_root_edge() const
{
  for (const auto& [eid, einfo] : graph.edges())
  {
    if (graph.incoming_edges(einfo.source).empty())
      return eid;
  }
  return -1;
}

}


#include "arbor.h"

namespace arbor::qsm {

QSM qsm(const PointCloud& tree, const settings::ArborParameters& params, const Logger& logger)
{
  size_t n = tree.size();

  // Filter wood only
  std::vector<bool> wood_mask(n, false);
  for (std::size_t i = 0; i < n; ++i) wood_mask[i] = tree.is_wood(i);
  PointCloud wood = tree.subset(wood_mask);

  // Determine the geographic coordinates minimum
  // and center on 0,0,0 for numerical stability
  size_t min_idx = 0;
  double min_z = std::numeric_limits<double>::max();
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

  // Sometime the very bottom have like two clusters of wood
  // this creates trouble. One of the two is remove to ensure
  // a single entry point for the QSM
  wood = QSM::clean_tree_butt(wood);

  // Inspired by aRchi and needed to build the skeleton
  auto layers = QSM::layers(wood, 0.2); //step = 0.2
  auto clusters = QSM::clusters(wood, layers, 0.1); // cl_dist = 0.1

  // Convert std::vector<pair>
  std::vector<std::pair<int, int>> iter_cluster;
  n = wood.size();
  iter_cluster.reserve(n);
  for (std::size_t i = 0; i < n; ++i) { iter_cluster.emplace_back(layers[i].first, clusters[i].first); }

  // Build the QSM
  QSM qsm;
  qsm.build_skeleton(wood, iter_cluster, 0.1); //max_d = 0.1
  qsm.compute_topology();

  // Fix duplicated root (rare)
  int n_root = 0;
  for (const auto& kv : qsm.cylinders()) {
    if (kv.second.parent_ID == 0)
      ++n_root;
  }

  if (n_root == 0) throw std::runtime_error("Internal error: 0 root");
  if (n_root > 1) qsm.fix_multiple_root();

  qsm.compute_architecture(1, false);
  qsm.smooth_skeleton(1, 0);
  qsm.detect_weird_butt();
  qsm.estimate_prolongation(wood);
  qsm.prolongate(qsm.prolongation_distance);
  qsm.construct_radii(wood);
  qsm.shift(tx, ty, tz);

  return qsm;
}

}

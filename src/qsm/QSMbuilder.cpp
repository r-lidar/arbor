#include "QSMbuilder.h"

void QSMbuilder::build(const PointCloud& tree)
{
  size_t n;

  // Filter wood only
  n = tree.size();
  std::vector<bool> wood_mask(n, false);
  for (std::size_t i = 0; i < n; ++i) wood_mask[i] = tree.is_wood(i);
  PointCloud wood = tree.subset(wood_mask);
  n = wood.size();

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

  // Sometime the very bottom have two clusters of wood
  // this creates trouble. One of the two is removed to ensure
  // a single entry point for the QSM
  wood = clean_tree_butt(wood);

  // Inspired by aRchi and needed to build the skeleton
  auto layers = this->layers(wood, 0.2); //step = 0.2
  auto clusters = this->clusters(wood, layers, 0.1); // cl_dist = 0.1

  // Convert std::vector<pair>
  std::vector<std::pair<int, int>> iter_cluster;
  iter_cluster.reserve(n);
  for (std::size_t i = 0; i < n; ++i) { iter_cluster.emplace_back(layers[i].first, clusters[i].first); }

  // Build the QSM nodes
  build_skeleton(wood, iter_cluster, 0.1); //max_d = 0.1

  // Connect the QSM nodes
  compute_topology();

  // Fix root issue (rare)
  int n_root = count_root();
  if (n_root == 0) throw std::runtime_error("Internal error: 0 root for this QSM. Please report.");
  if (n_root > 1) fix_multiple_root();

  compute_architecture(1, false);
  smooth_skeleton(1, 0);
  detect_weird_butt();
  estimate_prolongation(wood);
  prolongate(prolongation_distance);
  construct_radii(wood);
  shift(tx, ty, tz);
}

void QSMbuilder::shift(double tx, double ty, double tz)
{
  for (auto& kv : qsm)
  {
    kv.second.startX += tx;
    kv.second.endX += tx;
    kv.second.startY += ty;
    kv.second.endY += ty;
    kv.second.startZ += tz;
    kv.second.endZ += tz;
  }
}

int QSMbuilder::count_root()
{
  int n_root = 0;
  for (const auto& kv : qsm)
  {
    if (kv.second.parent_ID == 0)
      ++n_root;
  }
  return n_root;
}

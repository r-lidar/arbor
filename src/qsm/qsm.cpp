#include "arbor.h"

namespace arbor::qsm {

QSM qsm(const PointCloud& pc, const settings::ArborParameters& params, const Logger& logger)
{
  QSM qsm;

  auto layers = QSM::layers(pc, 0.2);
  auto clusters = QSM::clusters(pc, layers, 0.1);

  std::vector<std::pair<int, int>> iter_cluster;
  const std::size_t n = pc.size();
  iter_cluster.reserve(n);
  for (std::size_t i = 0; i < n; ++i) { iter_cluster.emplace_back(layers[i].first, clusters[i].first); }

  qsm.build_skeleton(pc, iter_cluster, 0.1);
  qsm.compute_topology();
  qsm.compute_architecture();
  qsm.smooth_skeleton();
  //qsm.detect_weird_butt(qsm);
  //qsm.estimate_prolongation(tree, qsm)
  //qsm.conic_allometry(qsm, 2*R0, tip_radius)
  //qsm.measure(tree, qsm, sarc = 180, sins = 0.2, sinl = 0.3, srmeas = 0.03)
  //qsm.polynomial_fitting(qsm, tip_radius)
  //qsm.reconstruction(qsm, tip_radius)
  return qsm;
}

}

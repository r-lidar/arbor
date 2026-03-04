#ifndef QSMBUILDER_H
#define QSMBUILDER_H

#include "arbor.h"
#include "QSM.h"

class QSMbuilder
{
public:
  QSMbuilder(QSM& qsm, const arbor::settings::ArborParameters& p = arbor::settings::ArborParameters()) : params(p), qsm(qsm) {};
  void build(const PointCloud& pc);
  void set_logger(Logger new_logger) { logger = std::move(new_logger); }

  // Static and public to be exported in R
  static std::vector<std::pair<int, double>> layers(const PointCloud& points, double D, const Logger& logger = [](const std::string&) {});
  static std::vector<std::pair<int, double>> clusters(const PointCloud& points, const std::vector<std::pair<int, double>>&, double cl_dist, const Logger& logger = [](const std::string&) {});
  static PointCloud clean_tree_butt(const PointCloud&, const Logger& logger = [](const std::string&) {});

  void build_skeleton(const PointCloud&, const std::vector<std::pair<int, int>>& iter_cluster, double max_d);
  void compute_topology();
  void compute_architecture(int root_id = 1, bool use_volume = true);
  void smooth_skeleton(int niter, double th);
  void detect_weird_butt(double thresh = 50.0, int window = 4);
  void estimate_prolongation(const PointCloud& tree);
  void prolongate(double d, double L = 0.1);
  void construct_radii(const PointCloud& tree, double tip_radius = 0.0025);
  void measure_radii(const PointCloud& tree, float sarc = 180, float sins = 0.2, float sinl = 0.3, float srmeas = 0.05);
  void polynomial_fitting(double tip_radius = 0.0025);
  void reconstruct_missing_radii(double tip_radius);
  void conic_allometry(double R0, double tip_radius = 0.0025);
  void fix_multiple_root();
  void shift(double tx, double ty, double tz);

  // recursive helpers
  double compute_subtree_length(int node_id);
  double compute_subtree_max_z(int node_id);
  double compute_subtree_volume(int node_id);
  void assign_subtree_ids(int node_id, int current_axis_id, int current_branch_order, int &next_axis_id, bool use_volume);

  // misc
  int count_root();
  void remove_disconnected_branches();
  double conic_allometry(double tip_radius, double wi, double w0, double r0) const;

  double prolongation_distance = 0;

  arbor::settings::ArborParameters params;
  QSM& qsm;

  Logger logger;
};

#endif

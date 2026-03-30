#ifndef QSMBUILDER_H
#define QSMBUILDER_H

#include "arbor.h"
#include "QSM.h"

namespace arbor::qsm {

class QSMbuilder
{
public:
  QSMbuilder(QSM& graph, const arbor::settings::ArborParameters& p = arbor::settings::ArborParameters()) : params(p), graph(graph) {};
  void build(const PointCloud& pc);

  // Static and public to be exported in R
  static std::vector<std::pair<int, double>> layers(const PointCloud& points, double D);
  static std::vector<std::pair<int, double>> clusters(const PointCloud& points, const std::vector<std::pair<int, double>>&, double cl_dist);
  static PointCloud clean_tree_butt(const PointCloud&);

  void build_skeleton(const PointCloud&, const std::vector<std::pair<int, int>>& iter_cluster, double max_d);
  void compute_topology();
  void compute_architecture(bool use_volume = false);
  void smooth_radii(int steps = 10, double lambda = 0.5, double mu = -0.7);
  void smooth_skeleton(int steps = 10, double lambda = 0.5, double mu = -0.53);
  void detect_weird_butt(double thresh = 50.0, int window = 4);
  void estimate_prolongation(const PointCloud& tree);
  void prolongate(double d, double L = 0.1);
  void construct_radii(const PointCloud& tree, double tip_radius = 0.0025);
  void measure_radii(const PointCloud& tree, float sarc = 180, float sins = 0.2, float sinl = 0.3, float srmeas = 0.05);
  void refine_radii(const PointCloud& tree);
  void polynomial_fitting(double tip_radius = 0.0025);
  void reconstruct_missing_radii(double tip_radius);
  void conic_allometry(double R0, double tip_radius = 0.0025);
  void fix_multiple_root();
  void shift(double tx, double ty, double tz);

  // recursive helpers (operate on graph edge IDs, analogous to cyl_IDs)
  double compute_subtree_length(int edge_id);
  double compute_subtree_max_z(int edge_id);
  double compute_subtree_volume(int edge_id);
  void assign_subtree_ids(int edge_id, int current_axis_id, int current_branch_order, int& next_axis_id, bool use_volume);

  // misc
  int count_nodes_connected_to_root() const;
  void remove_disconnected_branches();
  double conic_allometry(double tip_radius, double wi, double w0, double r0) const;

  // Find the root edge (first edge whose source node has no incoming edges).
  // Returns -1 if the graph is empty.
  int find_root_edge() const;

  double prolongation_distance = 0;

  arbor::settings::ArborParameters params;
  QSM& graph;
};

}

#endif

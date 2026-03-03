#ifndef QSM_H
#define QSM_H

#include "arbor.h"
#include "QSM.h"

class QSMbuilder
{
public:
  QSMbuilder(const arbor::settings::ArborParameter& p) : params(p)
  build(const PointCloud& pc);

private:
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
  void shift(double tx, double ty, double tz)
  {
    for (auto& kv : cylinders_)
    {
      kv.second.startX += tx;
      kv.second.endX += tx;
      kv.second.startY += ty;
      kv.second.endY += ty;
      kv.second.startZ += tz;
      kv.second.endZ += tz;
    }
  }

  // recursive helpers
  double compute_subtree_length(int node_id);
  double compute_subtree_max_z(int node_id);
  double compute_subtree_volume(int node_id);
  void assign_subtree_ids(int node_id, int current_axis_id, int current_branch_order, int &next_axis_id, bool use_volume);

private:
  arbor::settings::ArborParameter params;
  QSM qsm;
};

#endif

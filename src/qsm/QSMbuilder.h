/**
 * @file QSMbuilder.h
 * Project: Arbor
 *
 * Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

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

  // Static and public to be exposed in R
  static std::vector<int> cluster(const PointCloud& data, const std::vector<int>& iter, float eps);
  static std::unordered_map<int, std::vector<size_t>> group_points_by_edge(const QSM& graph, const PointCloud& tree);
  static PointCloud clean_tree_butt(const PointCloud&);

  void build_skeleton(const PointCloud&, const std::vector<std::pair<int, int>>& iter_cluster, double max_d);
  void compute_topology();
  void compute_architecture(bool use_volume = false);
  void smooth_radii();
  void smooth_skeleton(int steps = 10);
  void detect_weird_butt(double thresh = 50.0, int window = 4);
  void estimate_prolongation(const PointCloud& tree);
  void prolongate(double d, double L = 0.1);
  void construct_radii(const PointCloud& tree, double tip_radius = 0.0025);
  void measure_radii(const PointCloud& tree);
  void refine_radii(const PointCloud& tree);
  void refine_radii_broken(const PointCloud& tree);
  void polynomial_fitting(double tip_radius = 0.0025);
  void reconstruct_missing_radii(double tip_radius);
  void conic_allometry(double R0, double tip_radius = 0.0025);
  void fix_multiple_root();
  void shift(double tx, double ty, double tz);
  void prune_spurious_branches();
  void distance_to_root();

  // recursive helpers (operate on graph edge IDs, analogous to cyl_IDs)
  double compute_subtree_length(int edge_id);
  double compute_subtree_max_z(int edge_id);
  double compute_subtree_volume(int edge_id);
  void assign_subtree_ids(int edge_id, int current_axis_id, int current_branch_order, int& next_axis_id, bool use_volume);

  // misc
  int count_nodes_connected_to_root() const;
  void remove_disconnected_branches();
  double conic_allometry(double tip_radius, double wi, double w0, double r0) const;
  std::map<int, std::vector<int>> build_axis_map();

  // Find the root edge (first edge whose source node has no incoming edges).
  // Returns -1 if the graph is empty.
  int find_root_edge() const;

  double prolongation_distance = 0;

  arbor::settings::ArborParameters params;
  QSM& graph;
  bool likely_broken = false;
};

}

#endif

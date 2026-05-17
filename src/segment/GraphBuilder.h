/**
 * @file GraphBuilder.h
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

#ifndef GRAPHBUILDER_H
#define GRAPHBUILDER_H

#include "PointCloud.h"
#include "Graph.h"
#include "arbor.h"
#include "nanoflann.h"

#include <memory>

namespace arbor::segment {

class GraphBuilder
{
public:
  GraphBuilder(const settings::GraphParameters& p);
  ~GraphBuilder();
  Graph* get_graph();

  void add_core_layer(const PointCloud& core);
  void add_target_layer(const PointCloud& core, const PointCloud& target);
  void add_seed_layer(const PointCloud& core, const PointCloud& seeds);
  void add_master_seed_layer();
  void fix_directed_reachability(const PointCloud& cloud);

  void set_wood(const std::vector<bool>& x);

  int get_num_cores() const;
  int get_num_targets() const;
  int get_num_seeds() const;
  int get_num_master() const;
  std::pair<int, int> get_range_core() const;
  std::pair<int, int> get_range_targets() const;
  std::pair<int, int> get_range_seed() const;
  std::pair<int, int> get_range_master() const;

private:
  // KD-tree built once over the core point cloud in add_core_layer and shared
  // by add_target_layer, add_seed_layer, and fix_directed_reachability.
  // Released at the end of fix_directed_reachability.
  using KDTreeType = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 3>;

  Graph* graph;
  std::unique_ptr<KDTreeType> core_index;

  int offset_points  = -1;
  int offset_targets = -1;
  int offset_seeds   = -1;
  int offset_master  = -1;
  int total_core_nodes   = 0;
  int total_target_nodes = 0;
  int total_seed_nodes   = 0;
  int total_master_nodes = 0;
  int total_nodes        = 0;
  std::vector<bool> wood;
  bool graph_owner = true;
  settings::GraphParameters params;

  // Validates the angle penalty vector size; throws if invalid.
  void set_angle_penalty(const std::vector<float>& x);
};

}

#endif

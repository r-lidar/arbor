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

#include <memory>

namespace arbor::segment {

class GraphBuilder
{
public:
  GraphBuilder(const settings::GraphParameters& p);
  ~GraphBuilder();
  std::unique_ptr<Graph> get_graph();

  // Core kNN graph. Must be called first.
  void add_core_layer(const PointCloud& core);

  // Anchors each target on its nearest core point (Graph::target_anchors)
  void add_target_layer(const PointCloud& core, const PointCloud& target);

  // Connects each seed to its k_seed nearest core points (Graph::sources). For each core point
  // only the cheapest seed is kept, the first one in seed order on ties.
  void add_seed_layer(const PointCloud& core, const PointCloud& seeds);

  void fix_directed_reachability(const PointCloud& cloud);
  void fix_directed_reachability2(const PointCloud& cloud);
  void set_wood(const std::vector<bool>& x);

  int get_num_cores() const;
  int get_num_targets() const;
  int get_num_seeds() const;

private:
  std::unique_ptr<Graph> graph;

  int total_core_nodes   = 0;
  int total_target_nodes = 0;
  int total_seed_nodes   = 0;
  std::vector<bool> wood;
  settings::GraphParameters params;

  // Validates the angle penalty vector size; throws if invalid.
  void set_angle_penalty(const std::vector<float>& x);

  // Marks the core nodes reachable from the sources following directed edges
  std::vector<bool> reachable_from_sources() const;
};

}

#endif

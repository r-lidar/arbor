/**
 * @file segment_overfitting.h
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

#pragma once

#include <cstddef>
#include <vector>

#include "arbor.h"

namespace arbor::segment {

struct CircleRecord
{
  double x, y, z;
  double r;
  size_t n;         // number of distinct treeIDs supporting the circle
};

struct MatchRecord
{
  int ref;
  int move;
  size_t n;               // number of circles supporting this (ref, move) pair
  size_t ref_points  = 0; // cumulative inlier points attributed to `ref` across those circles
  size_t move_points = 0; // cumulative inlier points attributed to `move` across those circles

  // Ratio of the minority ID's total inlier points to the majority ID's,
  // aggregated across every circle that supports this match. Close to 1
  // means both IDs contributed comparable evidence; close to 0 means one
  // ID barely showed up (e.g. a couple of stray points from a neighboring
  // tree) and the match is likely a false positive rather than a genuine
  // over-segmented instance.
  double weight_ratio() const
  {
    if (ref_points == 0 || move_points == 0) return 0.0;
    const size_t lo = std::min(ref_points, move_points);
    const size_t hi = std::max(ref_points, move_points);
    return static_cast<double>(lo) / static_cast<double>(hi);
  }
};

struct InstanceMatchResult
{
  std::vector<CircleRecord> circles;
  std::vector<MatchRecord>  matches;
};

struct MergerConfig
{
  double layer_start      = 1.0;    // Min hag to make layers
  double layer_end        = 8.0;    // Max hag to make layers
  double layer_step       = 0.15;   // Layer thinkness
  float  dbscan_eps       = 0.05f;  // dbscan eps
  int    dbscan_min_pts   = 200;    // Min num points for dbscan cluster (speed)
  double fit_tolerance    = 0.02;   // Orbicular fitting tolerance (small)
  double min_arc_degree   = 240.0;  // Orbicular fitting minimal arc coverage for validity check
  double min_inlier_pct   = 60.0;   // Orbicular fitting minimal inlier percentage for validity check
  double max_interior_pct = 30.0;   // Orbicular fitting max insider percentage for validity check
  double max_radius       = 2.0;    // Orbicular fitting max radius (more than that it is not a tree)
  double min_radius       = 0.05;   // Orbicular fitting min radius (less than that there is not issue)
  int    min_support      = 3;      // Min number of OF fits to support a merging
  double min_weight_ratio = 0.1;    // Min percentage in the matching table.
};

// Detects trees that were split into multiple instances (over-segmentation)
// by scanning, layer by layer, for stem cross-sections whose inliers span
// more than one treeID, then merges those treeIDs back together.
class OverSegmentationResolver
{
public:
  explicit OverSegmentationResolver(MergerConfig config = MergerConfig());

  // Scans the point cloud and builds the circle/match evidence. Does not
  // modify `cloud`.
  InstanceMatchResult detect(const PointCloud& cloud) const;

  // Applies a matching table to `cloud`: for every (ref, move) pair supported
  // by more than config.min_support circles, remaps `move`'s treeID onto
  // `ref`. Chains of matches (e.g. 1-2, 2-3) are resolved via union-find
  // first, so every id in a chain collapses onto a single canonical id
  // instead of some remaps silently missing their target.
  void merge(PointCloud& cloud, const std::vector<MatchRecord>& matches) const;

  // Convenience: detect() followed by merge(). Returns the detection result
  // that was used to perform the merge.
  InstanceMatchResult run(PointCloud& cloud) const;

private:
  MergerConfig config_;
};

} // namespace arbor::segment

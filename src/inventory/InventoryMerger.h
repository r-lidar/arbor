/**
 * @file InventoryMerger.h
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

#ifndef INVENTORY_MERGER_H
#define INVENTORY_MERGER_H

#include "Inventory.h"
#include "PointCloud.h"
#include "fitting_types.h"

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

// Merging of a spatially accurate ground inventory into an instance segmentation
// (port of the arborInventoryMerger R package). The inventory is trusted: each
// inventory tree is searched where the inventory says it is, the stem is fitted
// at several heights around the point of measurement, and the IDs of the points
// that support the stem are mapped to the inventory ID.
//
// Pipeline (InventoryMerger::run() runs it all):
//   1. disentangle()  point cloud IDs that collide with inventory IDs are renumbered
//   2. match()        stem fitting for each inventory tree -> matching table (treeID, targetID)
//   3. reassign()     point cloud IDs are replaced by inventory IDs
//   4. resolve()      point cloud trees matched to several inventory trees are re-segmented
//
// Vocabulary: "treeID" or "point cloud ID" is the ID of arbor's segmentation, "targetID"
// or "inventory ID" is the ID of the inventory. Point cloud tree IDs <= 0 (including R's
// NA) mean "no tree" and are never remapped.

namespace arbor::inventory {

// Internal constants of the algorithm (hard-coded numbers in the R code)
struct MergerParameters
{
  std::vector<double> slice_offsets = {-2, -1.5, -1, -0.75, -0.5, 0, 0.5, 1, 1.5, 2}; // relative to the POM
  double slice_half_thickness   = 0.15;  // slice is [z+offset-h, z+offset+h]
  double min_slice_height       = 1.0;   // slices lower than that above ground are skipped
  double fit_tolerance          = 0.03;
  double min_arc_coverage       = 180;   // degrees. A slice fit is valid above this value
  double search_factor          = 1.5;   // initial search radius (x radius), incremented by the same amount on failure
  int    search_iterations      = 4;
  double neighbourhood_factor   = 5;     // neighbourhood of a tree is radius x this factor
  double neighbourhood_above    = 3;     // and up to z + this value
  double max_hag_above_pom      = 2.5;   // points above max(POM) + this value are never used
  double cluster_eps            = 0.1;   // dbscan of the neighbourhood
  int    cluster_min_pts        = 5;
  double cluster_max_base_hag   = 0.75;  // clusters that do not start below this height are discarded
};

// A stem cross-section successfully fitted for an inventory tree
struct StemSlice
{
  double offset = 0;                           // slice offset relative to the POM
  utils::fitting::FittingResult fit;           // center, radius, contour; inliers index 'points'
  std::vector<utils::fitting::Vec3> points;    // points of the slice
  std::vector<int> tree_ids;                   // unique point cloud IDs (> 0) of the inliers
};

// Everything found for one inventory tree.
// SingleCandidate: no stem fitted but a single point cloud tree stands at the inventory position.
struct StemMatch
{
  enum class Status { Fitted, SingleCandidate, NotDetected };

  int target_id = 0;
  Status status = Status::NotDetected;
  std::vector<StemSlice> slices;               // empty unless Fitted
  std::vector<int> tree_ids;                   // unique point cloud IDs to be merged into target_id
  std::vector<utils::fitting::Vec3> seeds;     // inliers of all slices or the inventory position

  int nfit() const { return tree_ids.empty() ? 0 : (status == Status::Fitted ? static_cast<int>(slices.size()) : 1); }
};

// One row of the matching table: the point cloud tree 'tree_id' is (part of) the inventory
// tree 'target_id'. tree_id = 0 means that no point cloud tree was matched (NA in R).
struct Match
{
  int tree_id = 0;
  int target_id = 0;
  int nfit = 0;
};

using MatchingTable = std::vector<Match>;

// One row (treeID, targetID, nfit) per point cloud tree matched, or (0, targetID, 0) if none
MatchingTable matching_table(const std::vector<StemMatch>& matches);

// Point cloud IDs matched to several distinct inventory IDs: treeID -> inventory IDs
std::map<int, std::vector<int>> find_conflicts(const MatchingTable& table);

struct MergeReport
{
  std::unordered_map<int, int> disentangled;   // point cloud ID renumbered in step 1: old -> new
  std::vector<StemMatch> matches;

  MatchingTable table() const { return matching_table(matches); }
  std::vector<int> undetected() const;         // inventory IDs with no point cloud tree
  std::map<int, std::vector<int>> conflicts() const { return find_conflicts(table()); }
  std::string to_string() const;
};

class InventoryMerger
{
public:
  // The inventory is not copied and must outlive the merger
  explicit InventoryMerger(const Inventory& inventory, MergerParameters params = {});

  // Full pipeline. Modifies the tree IDs of the point cloud in place and logs the report.
  MergeReport run(PointCloud& pc) const;

  // Renumbers the point cloud tree IDs that also exist in the inventory. Returns the map old -> new.
  std::unordered_map<int, int> disentangle(PointCloud& pc) const;

  // Fits the stems of all the inventory trees (in parallel). Throws if IDs are entangled.
  std::vector<StemMatch> match(const PointCloud& pc) const;

  // Fits the stems of some inventory trees only, e.g. to inspect or refit a tree. No ID check.
  std::vector<StemMatch> match(const PointCloud& pc, const std::vector<int>& ids) const;

  // Replaces point cloud IDs by inventory IDs. Point cloud IDs matched several times are left
  // untouched (see resolve()). Pre-existing IDs that collide with a target are moved to new IDs.
  static void reassign(PointCloud& pc, const MatchingTable& table);

  // Point cloud trees matched to several inventory trees are re-segmented with segment_instance()
  // using the stem inliers of each inventory tree as seeds. Requires 'hag' and 'foliage'.
  // 'matches' are those of match(); the inventory trees involved that are missing are fitted again.
  // Points not reached by the segmentation keep their ID.
  void resolve(PointCloud& pc, const MatchingTable& table, const std::vector<StemMatch>& matches = {}) const;

private:
  const Inventory& inventory_;
  MergerParameters params_;

  std::vector<StemMatch> fit(const PointCloud& pc, const std::vector<const InventoryTree*>& trees) const;
};

} // namespace arbor::inventory

#endif // INVENTORY_MERGER_H

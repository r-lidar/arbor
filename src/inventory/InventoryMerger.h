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

#include <string>
#include <unordered_map>
#include <vector>

namespace arbor::settings { struct ArborParameters; }

// Merging of a spatially accurate ground inventory into an instance segmentation
// (port of the arborInventoryMerger R package). The inventory is trusted: each
// inventory tree is searched where the inventory says it is, the stem is fitted
// at several heights around the point of measurement, and the IDs of the points
// that support the stem are mapped to the inventory ID.
//
// Pipeline (merge_inventory() runs it all):
//   1. disentangle_ids()        point cloud IDs that collide with inventory IDs are renumbered
//   2. matching_table()         stem fitting for each inventory tree -> (treeID, targetID) pairs
//   3. reassign_ids()           point cloud IDs are replaced by inventory IDs
//   4. resolve_multi_matching() point cloud trees matched to several inventory trees are re-segmented
//
// Point cloud tree IDs <= 0 (including R's NA) mean "no tree" and are never remapped.

namespace arbor::inventory {

// Internal constants of the algorithm. In the R code they are hard-coded numbers; the
// defaults are those numbers. They are grouped here only to give them a name.
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

// A stem cross-section successfully fitted for an inventory tree.
// R equivalent: one 'aim_fit' object returned by aim_fit_slice() (fit_circloid_cpp() output
// + to_be_merge_id + seeds).
struct SliceFit
{
  double offset = 0;                           // slice offset relative to the POM
  utils::fitting::FittingResult fit;           // center, radius, contour; inliers index 'points'
  std::vector<utils::fitting::Vec3> points;    // points of the slice
  std::vector<int> tree_ids;                   // unique point cloud IDs (> 0) of the inliers
};

// Everything found for one inventory tree.
// R equivalent: the list returned by aim_fit_tree(). 'tree_ids' is the union of the
// to_be_merge_id of all slices and 'seeds' the union of their seeds.
// Status::SingleCandidate and Status::NotDetected are the two warnings "No tree detected
// for tree X (but a single candidate ID was found and used)" of the R code.
struct TreeFit
{
  enum class Status { Fitted, SingleCandidate, NotDetected };

  int target_id = 0;
  Status status = Status::NotDetected;
  std::vector<SliceFit> slices;                // empty unless Fitted
  std::vector<int> tree_ids;                   // unique point cloud IDs to be merged into target_id
  std::vector<utils::fitting::Vec3> seeds;     // inliers of all slices or the inventory position

  int nfit() const { return tree_ids.empty() ? 0 : (status == Status::Fitted ? static_cast<int>(slices.size()) : 1); }
};

// One row of the matching table: the point cloud tree 'tree_id' is (part of) the inventory
// tree 'target_id'. R equivalent: a row of the data.frame returned by aim_matching_table().
// tree_id = 0 means that no point cloud tree was matched (NA in R).
struct Match
{
  int tree_id = 0;
  int target_id = 0;
  int nfit = 0;
};

using MatchingTable = std::vector<Match>;

struct MergeReport
{
  std::unordered_map<int, int> disentangled;   // point cloud ID renumbered in step 1: old -> new
  MatchingTable matching_table;
  std::vector<int> undetected;                 // inventory IDs with no point cloud tree
  std::unordered_map<int, std::vector<int>> conflicts; // point cloud ID -> several inventory IDs, re-segmented

  std::string to_string() const;
};

// Renumbers the point cloud tree IDs that also exist in the inventory. Returns the map old -> new.
std::unordered_map<int, int> disentangle_ids(PointCloud& pc, const Inventory& inventory);

// Fits the stem of one inventory tree. 'max_hag' bounds the points used (see MergerParameters::max_hag_above_pom).
TreeFit fit_tree(const PointCloud& pc, const InventoryTree& tree, double max_hag, const MergerParameters& params = {});

// Fits all the trees of the inventory (in parallel)
std::vector<TreeFit> fit_trees(const PointCloud& pc, const Inventory& inventory, const MergerParameters& params = {});

// Matching table from the fits, or fitting the trees first. Throws if IDs are entangled (see disentangle_ids())
MatchingTable matching_table(const std::vector<TreeFit>& fits);
MatchingTable matching_table(const PointCloud& pc, const Inventory& inventory, const MergerParameters& params = {});

// Replaces point cloud IDs by inventory IDs. Point cloud IDs matched several times are left
// untouched (see resolve_multi_matching()). Pre-existing IDs that collide with a target are moved
// to new IDs. Returns the number of point cloud IDs matched several times.
size_t reassign_ids(PointCloud& pc, const MatchingTable& table);

// Point cloud trees matched to several inventory trees are re-segmented with segment_instance()
// using the stem inliers of each inventory tree as seeds. Requires 'hag' and 'foliage'.
// Points not reached by the segmentation keep their ID. The version with fits reuses the fits
// of the matching step instead of re-fitting the trees.
void resolve_multi_matching(PointCloud& pc, const MatchingTable& table, const Inventory& inventory, const MergerParameters& params = {});
void resolve_multi_matching(PointCloud& pc, const MatchingTable& table, const std::vector<TreeFit>& fits);
void resolve_multi_matching(PointCloud& pc, const MatchingTable& table, const std::vector<TreeFit>& fits, const settings::ArborParameters& params);

// Full pipeline. Modifies the tree IDs of the point cloud in place and logs the report.
MergeReport merge_inventory(PointCloud& pc, const Inventory& inventory, const MergerParameters& params = {});

const char* to_string(TreeFit::Status status) noexcept;

} // namespace arbor::inventory

#endif // INVENTORY_MERGER_H

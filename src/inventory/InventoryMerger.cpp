/**
 * @file InventoryMerger.cpp
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

#include "InventoryMerger.h"

#include "arbor.h"
#include "myomp.h"
#include "nanoflann.h"
#include "progress.h"
#include "fitting_orbicular.h"
#include "dbscan/dbscan.hpp"

#include <algorithm>
#include <atomic>
#include <exception>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <unordered_set>

// Port of the R package arborInventoryMerger (AIM). Correspondence with the R code:
//
//   R (AIM)                          C++ (this file)
//   -------------------------------  ---------------------------------------------------
//   aim_disentangle_ids()            disentangle_ids()
//   aim_fit_tree()                   fit_tree()  (+ neighbourhood() for its first part)
//   aim_fit_slice()                  fit_slice()
//   aim_matching_table()             fit_trees() + matching_table()
//   aim_make_matching_table()        matching_table(const std::vector<TreeFit>&)
//   aim_reassign_ids()               reassign_ids()
//   aim_resolve_multi_matching()     resolve_multi_matching()
//   (README workflow, steps 5 to 8)  merge_inventory()
//
// Differences with the R code are flagged with "R difference" in the comments.
// Vocabulary: "treeID" or "point cloud ID" is the ID of arbor's segmentation, "targetID"
// or "inventory ID" is the ID of the inventory (ID_Arbre). IDs <= 0 (R's NA arrives as a
// large negative number) mean "no tree" and are never remapped.

namespace arbor::inventory {

using utils::fitting::Vec3;

namespace {

// ---------------------------------------------------------------------------
// 2D spatial index over a subset of the point cloud
//
// The R code scans the whole point cloud for every inventory tree:
//   pc <- las@data[sqrt((X - xc)^2 + (Y - yc)^2) < radius*5 & ...]
// Here a KD-tree on X,Y is built once and each tree queries only its neighbourhood.
// Only the points that can be useful are indexed (hag <= max(POM) + 2.5, like the
// 'useful' point cloud of aim_matching_table()), which keeps the index small.
// ---------------------------------------------------------------------------

struct XYAdaptor
{
  const PointCloud& pc;
  std::vector<unsigned int> idx; // indices of the indexed points in pc

  inline size_t kdtree_get_point_count() const { return idx.size(); }
  inline double kdtree_get_pt(size_t i, size_t dim) const { return dim == 0 ? pc.get_x(idx[i]) : pc.get_y(idx[i]); }
  template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
};

using XYTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, XYAdaptor>, XYAdaptor, 2>;

struct SpatialIndex
{
  XYAdaptor adaptor;
  XYTree tree;

  // Indexes the points with hag <= max_hag
  SpatialIndex(const PointCloud& pc, double max_hag) : adaptor{pc, {}}, tree(2, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10))
  {
    for (size_t i = 0; i < pc.size(); ++i)
    {
      if (pc.get_hag(i) <= max_hag)
        adaptor.idx.push_back(static_cast<unsigned int>(i));
    }
    tree.buildIndex();
  }

  // Indices in pc of the points within a radius in XY
  std::vector<unsigned int> query(double x, double y, double radius) const
  {
    std::vector<nanoflann::ResultItem<XYTree::IndexType, double>> matches;
    double q[2] = {x, y};
    tree.radiusSearch(q, radius * radius, matches, nanoflann::SearchParameters(0.0f, false));

    std::vector<unsigned int> out;
    out.reserve(matches.size());
    for (const auto& m : matches) out.push_back(adaptor.idx[m.first]);
    std::sort(out.begin(), out.end());
    return out;
  }
};

// A point of the neighbourhood of an inventory tree. Equivalent of a row of the R
// data.table 'pc' with columns X, Y, Z, hag, treeID. Points are copied because the
// neighbourhood of a tree is small and it is accessed many times.
struct NPoint
{
  double x, y, z, hag;
  int id;
};

// false for "no tree" (0, -1 or R's NA)
inline bool is_tree(int id) { return id > 0; }

// Equivalent of R unique() keeping the order of first appearance

template <typename T>
void push_unique(std::vector<T>& v, T x)
{
  if (std::find(v.begin(), v.end(), x) == v.end()) v.push_back(x);
}

void check_pointcloud(const PointCloud& pc, const char* fun)
{
  if (!pc.has_treeid()) throw std::runtime_error(std::string(fun) + ": point cloud is missing required 'treeID' attribute.");
  if (!pc.has_hag())    throw std::runtime_error(std::string(fun) + ": point cloud is missing required 'hag' attribute.");
}

// R: z <- max(inventory$POM) + 2.5 ; useful <- filter_poi(las, hag <= z)
double max_hag_of(const Inventory& inventory, const MergerParameters& params)
{
  double max_pom = 0;
  for (const InventoryTree& t : inventory) max_pom = std::max(max_pom, t.pom);
  return max_pom + params.max_hag_above_pom;
}

// First part of R aim_fit_tree(): the points around the inventory tree.
//
//   pc <- las@data[sqrt((X - xc)^2 + (Y - yc)^2) < radius*5 & Z < zc+3]
//   cl <- dbscan::dbscan(pc[, .(X,Y,Z)], eps = 0.1)
//   min_height_by_cl = pc[cl > 0, min(hag), by = cl]
//   keep_cl = min_height_by_cl[V1 < 0.75]$cl
//   pc = pc[cl %in% keep_cl]
//
// i.e. the points within 5 radii of the inventory position and below 3 m above the POM,
// clustered with dbscan. Only the clusters that touch the ground (min hag < 0.75) are
// kept: this removes branches and crowns of the neighbours hanging around the stem.
// Noise points (R cluster 0) are not part of any cluster and are thus removed too.
std::vector<NPoint> neighbourhood(const PointCloud& pc, const SpatialIndex& index, const InventoryTree& tree, const MergerParameters& params)
{
  // R: las@data[sqrt((X - xc)^2 + (Y - yc)^2) < radius*5 & Z < zc+3]
  std::vector<NPoint> pts;
  for (unsigned int i : index.query(tree.x, tree.y, tree.radius * params.neighbourhood_factor))
  {
    double z = pc.get_z(i);
    if (z < tree.z + params.neighbourhood_above)
      pts.push_back({pc.get_x(i), pc.get_y(i), z, pc.get_hag(i), pc.get_treeid(i)});
  }

  if (pts.empty()) return pts;

  // Coordinates relative to the tree for numerical accuracy
  std::vector<point3> xyz;
  xyz.reserve(pts.size());
  for (const NPoint& p : pts) xyz.push_back({p.x - tree.x, p.y - tree.y, p.z - tree.z});

  // R: dbscan::dbscan(pc[, .(X,Y,Z)], eps = 0.1), minPts = 5 is the default of the R package.
  // R difference: the vendored dbscan may attribute a border point reachable from two
  // clusters differently. This has no practical consequence.
  std::vector<std::vector<size_t>> clusters = dbscan(xyz, static_cast<float>(params.cluster_eps), params.cluster_min_pts);

  // R: keep the clusters whose lowest point is below 0.75 m above ground
  std::vector<NPoint> kept;
  for (const std::vector<size_t>& cluster : clusters)
  {
    double min_hag = std::numeric_limits<double>::infinity();
    for (size_t k : cluster) min_hag = std::min(min_hag, pts[k].hag);
    if (min_hag >= params.cluster_max_base_hag) continue;
    for (size_t k : cluster) kept.push_back(pts[k]);
  }

  return kept;
}

// R aim_fit_slice(): fits the stem cross-section in a 30 cm thick horizontal slice
// centered at POM + offset. Returns false where R returns NULL (no valid stem section).
// The candidates are the points of the neighbourhood within the current search radius.
bool fit_slice(const std::vector<NPoint>& candidates, const InventoryTree& tree, double offset, const MergerParameters& params, SliceFit& out)
{
  // R: if (zc + zoffset < 1) return(NULL)
  // R difference (bug fix): zc is an absolute elevation, so the R test practically never
  // triggered. The intent is to skip slices lower than 1 m above ground, i.e. POM + offset.
  if (tree.pom + offset < params.min_slice_height) return false;

  double zc = tree.z + offset;
  double zmin = zc - params.slice_half_thickness;
  double zmax = zc + params.slice_half_thickness;

  // R: slice = slice[Z > zc + zoffset - 0.15 & Z < zc + zoffset + 0.15]
  std::vector<const NPoint*> slice;
  double min_hag = std::numeric_limits<double>::infinity();
  double max_hag = -std::numeric_limits<double>::infinity();
  for (const NPoint& p : candidates)
  {
    if (p.z > zmin && p.z < zmax)
    {
      slice.push_back(&p);
      min_hag = std::min(min_hag, p.hag);
      max_hag = std::max(max_hag, p.hag);
    }
  }

  // R: same three tests. Not enough points, slice entirely close to the ground, or slice
  // entirely much higher than the POM (the stem does not reach the ground: not our tree)
  if (slice.size() < 3) return false;
  if (max_hag < 0.5) return false;
  if (min_hag > 2 * tree.pom) return false;

  // R: res <- arbor:::fit_circloid_cpp(xyz, tolerance = 0.03, complexity = 3)
  // fit_circloid_cpp() with its default axis uses a vertical axis through the centroid of
  // the points, and complexity = 3 is FitMode::Buttress. This is reproduced exactly.
  double mx = 0, my = 0;
  for (const NPoint* p : slice) { mx += p->x; my += p->y; }
  mx /= slice.size();
  my /= slice.size();

  utils::fitting::CrossSectionFitter fitter;
  fitter.set_axis({mx, my, 0}, {mx, my, 1});
  for (const NPoint* p : slice) fitter.add_point(p->x, p->y, p->z);

  utils::fitting::FittingResult res = fitter.fit(params.fit_tolerance, utils::fitting::FitMode::Buttress);
  // R: if (res$covered_arc_degree > 180) keep the fit, otherwise NULL. A fit that covers
  // less than half of the circle is not trusted to be a stem.
  if (!res.success || res.arc_coverage_deg <= params.min_arc_coverage) return false;

  // R: res$to_be_merge_id = unique(inliner$treeID) ; res$seeds = slice[res$inliers, .(X,Y,Z)]
  // The points of the slice are kept so that the inliers (indices in 'points') can be
  // converted into seeds later. R difference: points without tree (NA) are not included in
  // the IDs to be merged, in R they produce an NA that is dropped later anyway.

  out.offset = offset;
  out.points.clear();
  out.points.reserve(slice.size());
  for (const NPoint* p : slice) out.points.emplace_back(p->x, p->y, p->z);
  out.tree_ids.clear();
  for (int k : res.inlier_indices)
  {
    if (is_tree(slice[k]->id)) push_unique(out.tree_ids, slice[k]->id);
  }
  out.fit = std::move(res);
  return true;
}

// R aim_fit_tree(): searches the stem of one inventory tree and returns the fits, the
// point cloud IDs to be merged into the inventory tree and the seeds (stem points) that
// resolve_multi_matching() uses to re-segment the trees.
//
// The warning is returned instead of being emitted because this function runs in OpenMP
// threads, from which messages are discarded (see progress.h).
TreeFit fit_tree(const PointCloud& pc, const SpatialIndex& index, const InventoryTree& tree, const MergerParameters& params, std::string& warning)
{
  TreeFit res;
  res.target_id = tree.id;

  std::vector<NPoint> pts = neighbourhood(pc, index, tree, params);

  // R:
  //   while (it <= 3) {
  //     slice = pc[sqrt((X - xc)^2 + (Y - yc)^2) < radius*search_radius_factor]
  //     for (i in seq_along(offsets)) res[[i]] <- aim_fit_slice(slice, ..., offsets[i])
  //     if (length(res) == 0) search_radius_factor = search_radius_factor + 1.5 else stop
  //   }
  // The stem is searched at 10 heights around the POM in a disc of 1.5 radius. If nothing
  // is found, the disc grows to 3, 4.5 and 6 radii (4 iterations in total).
  std::vector<NPoint> candidates;
  double factor = params.search_factor;
  for (int it = 0; it < params.search_iterations && res.slices.empty(); ++it, factor += params.search_factor)
  {
    double r2 = (tree.radius * factor) * (tree.radius * factor);
    candidates.clear();
    for (const NPoint& p : pts)
    {
      double dx = p.x - tree.x, dy = p.y - tree.y;
      if (dx*dx + dy*dy < r2) candidates.push_back(p);
    }

    for (double offset : params.slice_offsets)
    {
      SliceFit slice;
      if (fit_slice(candidates, tree, offset, params, slice))
        res.slices.push_back(std::move(slice));
    }
  }

  // At least one slice fitted: the IDs to merge are the union of the IDs of all slices
  // (R aim_make_matching_table(): unique(do.call(c, to_be_merge_id))) and the seeds are
  // the inliers of all slices (R aim_resolve_multi_matching(): do.call(rbind, seeds)).
  if (!res.slices.empty())
  {
    res.status = TreeFit::Status::Fitted;
    for (const SliceFit& slice : res.slices)
    {
      for (int id : slice.tree_ids) push_unique(res.tree_ids, id);
      for (int k : slice.fit.inlier_indices) res.seeds.push_back(slice.points[k]);
    }
    return res;
  }

  // No stem found. R:
  //   slice <- slice[Z > zc - 0.5 & Z < zc + 0.5 & (X - xc)^2 < 2*radius & (Y - yc)^2 < 2*radius]
  //   candidates <- unique(slice$treeID)
  //   if (length(candidates) == 1) use it (with a warning) else NA (with a warning)
  // If a single point cloud tree is present at the inventory position it is used anyway.
  // The seed is then the inventory position itself (R: seeds = xyz).
  // R difference (bug fix): the R test compares a squared distance with a length in a
  // square box. It is now a disc of 2 radii. Points without tree are not candidates.
  res.seeds.emplace_back(tree.x, tree.y, tree.z);

  std::vector<int> ids;
  double r2 = 4 * tree.radius * tree.radius;
  for (const NPoint& p : candidates)
  {
    double dx = p.x - tree.x, dy = p.y - tree.y;
    if (p.z > tree.z - 0.5 && p.z < tree.z + 0.5 && dx*dx + dy*dy < r2 && is_tree(p.id))
      push_unique(ids, p.id);
  }

  if (ids.size() == 1)
  {
    res.status = TreeFit::Status::SingleCandidate;
    res.tree_ids = ids;
    warning = "No tree detected for tree " + std::to_string(tree.id) + " but a single candidate ID was found and used.";
  }
  else
  {
    res.status = TreeFit::Status::NotDetected;
    warning = "No tree detected for tree " + std::to_string(tree.id);
  }

  return res;
}

// The loop of R aim_matching_table(): for (i in seq_len(n)) aim_fit_tree(useful, inventory[i,])
// Fits the trees at positions 'which' in the inventory. The trees are independent and are
// processed in parallel. Warnings are collected per tree and emitted after the loop, in
// inventory order. An exception in a thread stops the loop and is rethrown afterwards.
std::vector<TreeFit> fit_trees(const PointCloud& pc, const Inventory& inventory, const std::vector<size_t>& which, const MergerParameters& params)
{
  check_pointcloud(pc, "fit_trees");

  std::vector<TreeFit> fits(which.size());
  if (which.empty()) return fits;

  SpatialIndex index(pc, max_hag_of(inventory, params));

  std::vector<std::string> warnings(which.size());
  std::exception_ptr error = nullptr;
  std::atomic<bool> failed{false};

  Counter counter("Fitting inventory trees", which.size());

  #pragma omp parallel for schedule(dynamic)
  for (size_t k = 0; k < which.size(); ++k)
  {
    if (counter.cancelled() || failed.load(std::memory_order_relaxed)) continue;

    try
    {
      fits[k] = fit_tree(pc, index, inventory[which[k]], params, warnings[k]);
    }
    catch (...)
    {
      #pragma omp critical
      if (!error) error = std::current_exception();
      failed = true;
    }

    counter.tick();
  }

  counter.throw_if_cancelled();
  if (error) std::rethrow_exception(error);

  for (const std::string& w : warnings)
    if (!w.empty()) logging::warning(w);

  return fits;
}

// Point cloud IDs matched to several distinct inventory IDs: the "conflicts" printed by
// print.aim_matching_table(), the duplicated treeIDs of aim_reassign_ids() and the dupIDs
// of aim_resolve_multi_matching(). Typically one arbor tree that is in reality two or more
// inventory trees. Returns treeID -> inventory IDs.
// R difference: R counts duplicated rows, here identical (treeID, targetID) rows given
// twice are not a conflict.
std::map<int, std::vector<int>> find_conflicts(const MatchingTable& table)
{
  std::map<int, std::vector<int>> targets;
  for (const Match& m : table)
  {
    if (is_tree(m.tree_id)) push_unique(targets[m.tree_id], m.target_id);
  }

  std::map<int, std::vector<int>> conflicts;
  for (auto& [id, t] : targets)
  {
    if (t.size() > 1) conflicts.emplace(id, std::move(t));
  }
  return conflicts;
}

// R aim_reassign_ids(). The matching table gives rules treeID -> targetID. The point cloud
// is relabelled accordingly, taking care that no two trees end up merged by accident.
//
// The R code works in 3 phases on the vector of IDs:
//   Phase 0: pre-existing IDs that are equal to a target (and that are not themselves a
//            source) are moved to fresh IDs, otherwise they would merge with the points
//            relabelled to this target.
//   Phase 1: source -> temporary ID
//   Phase 2: temporary ID -> target
// The final ID of a point only depends on its original ID, so the 3 phases are collapsed
// into a single lookup table old ID -> new ID applied once. Fresh IDs are the same as in R.
//
// 'warn' is false in merge_inventory() where conflicts are expected and handled.
size_t reassign_ids(PointCloud& pc, const MatchingTable& table, bool warn)
{
  if (!pc.has_treeid()) throw std::runtime_error("reassign_ids: point cloud is missing required 'treeID' attribute.");

  for (const Match& m : table)
  {
    if (m.target_id <= 0) throw std::invalid_argument("reassign_ids: invalid target ID " + std::to_string(m.target_id) + " in the matching table.");
  }

  std::map<int, std::vector<int>> conflicts = find_conflicts(table);
  if (warn && !conflicts.empty())
    logging::warning("Matching table contains duplicated 'treeID': use resolve_multi_matching().");

  // R: rows with NA treeID removed, rows with duplicated treeID removed (they are handled by
  // resolve_multi_matching()), identity rows (54 -> 54) removed.
  // Rules source -> target. Identity rules change nothing and must not be considered as handled
  // IDs, otherwise they are not protected from the points reassigned to them.
  std::vector<int> src, tgt;
  std::unordered_set<int> seen;
  for (const Match& m : table)
  {
    if (!is_tree(m.tree_id) || m.tree_id == m.target_id) continue;
    if (conflicts.count(m.tree_id) || !seen.insert(m.tree_id).second) continue;
    src.push_back(m.tree_id);
    tgt.push_back(m.target_id);
  }

  if (src.empty()) return conflicts.size();

  // Targets that are not sources collide with the pre-existing IDs of the point cloud: these
  // IDs are moved to new IDs. A target that is also a source is freed when its source moves.
  std::unordered_set<int> sources(src.begin(), src.end());
  std::unordered_set<int> colliding;
  for (int t : tgt)
    if (!sources.count(t)) colliding.insert(t);

  // R: bump_values <- intersect(unique(v), colliding_targets)
  //    max_val <- max(c(v, src_ids, tgt_ids))
  int max_val = 0;
  std::vector<int> bump; // pre-existing IDs to move, in order of appearance
  std::unordered_set<int> bumped;
  for (size_t i = 0; i < pc.size(); ++i)
  {
    int id = pc.get_treeid(i);
    if (!is_tree(id)) continue;
    max_val = std::max(max_val, id);
    if (colliding.count(id) && bumped.insert(id).second) bump.push_back(id);
  }
  for (int id : src) max_val = std::max(max_val, id);
  for (int id : tgt) max_val = std::max(max_val, id);

  // R: new_pool <- max_val + seq_len(n_src + n_bump). The first n_src IDs are the temporary
  // IDs of phase 1 (not needed here), the next ones are the fresh IDs of phase 0.
  std::unordered_map<int, int> map;
  for (size_t k = 0; k < src.size(); ++k) map[src[k]] = tgt[k];
  for (size_t k = 0; k < bump.size(); ++k) map[bump[k]] = max_val + static_cast<int>(src.size() + k + 1);

  for (size_t i = 0; i < pc.size(); ++i)
  {
    auto it = map.find(pc.get_treeid(i));
    if (it != map.end()) pc.set_treeid(i, it->second);
  }

  return conflicts.size();
}

} // namespace

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

const char* to_string(TreeFit::Status status) noexcept
{
  switch (status)
  {
  case TreeFit::Status::Fitted:          return "fitted";
  case TreeFit::Status::SingleCandidate: return "single candidate";
  case TreeFit::Status::NotDetected:     return "not detected";
  }
  return "unknown";
}

// R aim_disentangle_ids(). If a point cloud tree has an ID that also exists in the
// inventory (e.g. arbor's tree 54 and inventory's tree 54, almost surely two different
// trees), the point cloud tree gets a new ID. After that, point cloud IDs and inventory
// IDs never overlap and reassign_ids() cannot confuse them.
std::unordered_map<int, int> disentangle_ids(PointCloud& pc, const Inventory& inventory)
{
  if (!pc.has_treeid()) throw std::runtime_error("disentangle_ids: point cloud is missing required 'treeID' attribute.");

  // Point cloud IDs, in order of appearance
  std::vector<int> ids;
  std::unordered_set<int> seen;
  int max_id = 0;
  for (size_t i = 0; i < pc.size(); ++i)
  {
    int id = pc.get_treeid(i);
    if (!is_tree(id)) continue;
    if (seen.insert(id).second) ids.push_back(id);
    max_id = std::max(max_id, id);
  }

  std::unordered_map<int, int> map;

  // R: new_ids <- seq(from = max(las_ids) + 1, ...)
  // R difference (bug fix): new IDs start after the max of both the point cloud and the
  // inventory IDs. In R a new ID could be equal to an inventory ID.
  int next = std::max(max_id, inventory.next_id() - 1) + 1;
  for (int id : ids)
  {
    if (inventory.contains(id)) map[id] = next++;
  }

  if (map.empty())
  {
    logging::info("No matching treeIDs found between the point cloud and the inventory.");
    return map;
  }

  for (size_t i = 0; i < pc.size(); ++i)
  {
    auto it = map.find(pc.get_treeid(i));
    if (it != map.end()) pc.set_treeid(i, it->second);
  }

  return map;
}

// R aim_fit_tree() as called by the user for debugging (arborInventoryMerger:::aim_fit_tree).
// Builds a spatial index for a single tree, which is fine for occasional calls.
TreeFit fit_tree(const PointCloud& pc, const InventoryTree& tree, double max_hag, const MergerParameters& params)
{
  check_pointcloud(pc, "fit_tree");
  SpatialIndex index(pc, max_hag);
  std::string warning;
  TreeFit res = fit_tree(pc, index, tree, params, warning);
  if (!warning.empty()) logging::warning(warning);
  return res;
}

std::vector<TreeFit> fit_trees(const PointCloud& pc, const Inventory& inventory, const MergerParameters& params)
{
  std::vector<size_t> all(inventory.size());
  for (size_t i = 0; i < all.size(); ++i) all[i] = i;
  return fit_trees(pc, inventory, all, params);
}

// R aim_make_matching_table() applied to each tree: one row (treeID, targetID, nfit) per
// point cloud tree to merge into the inventory tree, or a single row (NA, targetID, 0) if
// nothing was found. nfit is the number of slices successfully fitted (1 for the single
// candidate fallback).
MatchingTable matching_table(const std::vector<TreeFit>& fits)
{
  MatchingTable table;
  for (const TreeFit& fit : fits)
  {
    if (fit.tree_ids.empty())
    {
      table.push_back({0, fit.target_id, 0});
      continue;
    }

    for (int id : fit.tree_ids)
      table.push_back({id, fit.target_id, fit.nfit()});
  }
  return table;
}

// R aim_matching_table()
MatchingTable matching_table(const PointCloud& pc, const Inventory& inventory, const MergerParameters& params)
{
  check_pointcloud(pc, "matching_table");

  // R: if (any(treeID %in% targetID)) stop(...)
  for (size_t i = 0; i < pc.size(); ++i)
  {
    if (inventory.contains(pc.get_treeid(i)))
      throw std::runtime_error("Some tree IDs in the point cloud are also in the inventory tree IDs. This can cause matching issues. Call disentangle_ids() first.");
  }

  return matching_table(fit_trees(pc, inventory, params));
}

size_t reassign_ids(PointCloud& pc, const MatchingTable& table)
{
  return reassign_ids(pc, table, true);
}

// R aim_resolve_multi_matching(). When one arbor tree (e.g. 123) is matched to several
// inventory trees (e.g. 12 and 13), the points of tree 123 are re-segmented with arbor's
// segment_instance() using the stem points found for 12 and 13 as seeds. Each point then
// goes to the closest seed in the graph sense, which splits the tree.
//
// R:
//   seeds <- inliers of aim_fit_tree() for each inventory tree involved, treeID = ID_Arbre
//   dupTrees <- filter_poi(las, treeID %in% dupIDs)
//   dupTreesFixed <- arbor::segment_instance(dupTrees, seeds)
//   las$treeID[dupTreesFixed$PID] = dupTreesFixed$treeID
//
// All the conflicting trees are re-segmented together in a single segment_instance() call.
void resolve_multi_matching(PointCloud& pc, const MatchingTable& table, const std::vector<TreeFit>& fits, const settings::ArborParameters& params)
{
  std::map<int, std::vector<int>> conflicts = find_conflicts(table);
  if (conflicts.empty()) return;

  if (!pc.has_treeid())  throw std::runtime_error("resolve_multi_matching: point cloud is missing required 'treeID' attribute.");
  if (!pc.has_hag())     throw std::runtime_error("resolve_multi_matching: point cloud is missing required 'hag' attribute.");
  if (!pc.has_foliage()) throw std::runtime_error("resolve_multi_matching: point cloud is missing required 'foliage' attribute.");

  std::unordered_map<int, const TreeFit*> fit_of;
  for (const TreeFit& fit : fits) fit_of[fit.target_id] = &fit;

  // Seeds: stem inliers of each inventory tree involved. R difference: the fits come from the
  // matching step instead of being recomputed on the full point cloud (see the overload below)
  std::vector<std::pair<Vec3, int>> seed_pts;
  for (const auto& [id, targets] : conflicts)
  {
    for (int target : targets)
    {
      auto it = fit_of.find(target);
      if (it == fit_of.end()) throw std::runtime_error("resolve_multi_matching: no fit for inventory tree " + std::to_string(target));
      for (const Vec3& p : it->second->seeds) seed_pts.emplace_back(p, target);
    }
  }

  // R: dupTrees <- filter_poi(las, treeID %in% dupIDs). 'points' plays the role of the R PID
  // column: the index of each point of the subset in the original point cloud.
  std::vector<unsigned int> points;
  for (size_t i = 0; i < pc.size(); ++i)
  {
    if (conflicts.count(pc.get_treeid(i))) points.push_back(static_cast<unsigned int>(i));
  }

  if (points.empty() || seed_pts.empty()) return;

  // segment_instance() reorders the points (partition above/below cut_above_ground), so
  // the position of a point in the subset changes. The R code relies on the PID column,
  // which is moved with the points. Here the original index is stored in the 'passage'
  // attribute, which segment_instance() moves with the points and does not use.
  PointCloud sub(points.size(), true);
  for (size_t k = 0; k < points.size(); ++k)
  {
    unsigned int i = points[k];
    sub.set_x(k, pc.get_x(i));
    sub.set_y(k, pc.get_y(i));
    sub.set_z(k, pc.get_z(i));
    sub.set_hag(k, pc.get_hag(i));
    sub.set_foliage(k, pc.get_foliage(i));
    sub.set_treeid(k, 0);
    sub.set_passage(k, static_cast<int>(k));
  }

  PointCloud seeds(seed_pts.size(), true);
  for (size_t k = 0; k < seed_pts.size(); ++k)
  {
    seeds.set_x(k, seed_pts[k].first.x);
    seeds.set_y(k, seed_pts[k].first.y);
    seeds.set_z(k, seed_pts[k].first.z);
    seeds.set_treeid(k, seed_pts[k].second);
  }

  segment::segment_instance(sub, seeds, params);

  // R: las$treeID[dupTreesFixed$PID] = dupTreesFixed$treeID
  // R difference: points not reached by the segmentation (NA) or below cut_above_ground (0)
  // keep their original ID instead of losing it. true_size() includes the points moved
  // beyond size() by the partition.
  for (size_t k = 0; k < sub.true_size(); ++k)
  {
    int id = sub.get_treeid(k);
    if (is_tree(id)) pc.set_treeid(points[sub.get_passage(k)], id);
  }
}

// Same with default arbor parameters (the R code uses arbor_parameters_default)
// R difference: the oversegmentation solver is disabled.
void resolve_multi_matching(PointCloud& pc, const MatchingTable& table, const std::vector<TreeFit>& fits)
{
  settings::ArborParameters params;
  params.instance.oversegmentation_solver_enabled = false; // it would merge the trees we are splitting
  resolve_multi_matching(pc, table, fits, params);
}

// R aim_resolve_multi_matching(las, matching_table, inventory) as called from R: the fits of
// the matching step are not available, so only the inventory trees involved in a conflict
// are fitted again.
// R difference: the fits use the same point cloud filter as the matching step (hag <= max(POM)
// + 2.5) instead of the full point cloud, which is faster and consistent with the matching.
void resolve_multi_matching(PointCloud& pc, const MatchingTable& table, const Inventory& inventory, const MergerParameters& params)
{
  std::map<int, std::vector<int>> conflicts = find_conflicts(table);
  if (conflicts.empty()) return;

  std::vector<size_t> which;
  std::unordered_set<int> done;
  for (const auto& [id, targets] : conflicts)
  {
    for (int target : targets)
      if (done.insert(target).second) which.push_back(inventory.index_of(target));
  }

  std::vector<TreeFit> fits = fit_trees(pc, inventory, which, params);
  resolve_multi_matching(pc, table, fits);
}

// The whole README workflow of arborInventoryMerger, for use outside R (e.g. a Qt GUI):
//
//   las <- aim_disentangle_ids(las, inventory)
//   matching_table <- aim_matching_table(las, inventory)
//   las$treeID <- aim_reassign_ids(las, matching_table)
//   las$treeID <- aim_resolve_multi_matching(las, matching_table, inventory)
//
// The tree IDs are modified in place. The fits of the matching step are reused to resolve
// the multi-matching instead of being recomputed.
MergeReport merge_inventory(PointCloud& pc, const Inventory& inventory, const MergerParameters& params)
{
  check_pointcloud(pc, "merge_inventory");
  if (inventory.empty()) throw std::runtime_error("merge_inventory: inventory is empty.");

  MergeReport report;

  Task task("Inventory merging", {1, 20, 1, 10});

  task.step("Disentangling IDs");
  report.disentangled = disentangle_ids(pc, inventory);

  task.step("Matching inventory trees");
  std::vector<TreeFit> fits = fit_trees(pc, inventory, params);
  report.matching_table = matching_table(fits);

  for (const Match& m : report.matching_table)
    if (!is_tree(m.tree_id)) report.undetected.push_back(m.target_id);
  for (auto& [id, targets] : find_conflicts(report.matching_table))
    report.conflicts.emplace(id, std::move(targets));

  task.step("Reassigning IDs");
  reassign_ids(pc, report.matching_table, false);

  task.step("Resolving multi-matching");
  resolve_multi_matching(pc, report.matching_table, fits);

  logging::info(report.to_string());

  return report;
}

// Human readable summary, similar to print.aim_matching_table()
std::string MergeReport::to_string() const
{
  std::set<int> targets, matched;
  for (const Match& m : matching_table)
  {
    targets.insert(m.target_id);
    if (m.tree_id > 0) matched.insert(m.target_id);
  }

  std::ostringstream os;
  os << "Inventory merging: " << targets.size() << " inventory trees | "
     << matched.size() << " matched | "
     << undetected.size() << " undetected | "
     << conflicts.size() << " conflicts | "
     << disentangled.size() << " IDs renumbered";

  if (!undetected.empty())
  {
    os << "\nUndetected:";
    for (int id : undetected) os << " " << id;
  }

  if (!conflicts.empty())
  {
    std::map<int, std::vector<int>> sorted(conflicts.begin(), conflicts.end());
    os << "\nConflicts (re-segmented):";
    for (const auto& [id, t] : sorted)
    {
      os << " " << id << "->";
      for (size_t k = 0; k < t.size(); ++k) os << (k ? "," : "") << t[k];
    }
  }

  return os.str();
}

} // namespace arbor::inventory

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
#include <optional>
#include <set>
#include <sstream>
#include <unordered_set>

// Port of the R package arborInventoryMerger (AIM):
//   aim_disentangle_ids()         InventoryMerger::disentangle()
//   aim_matching_table()          InventoryMerger::match() + matching_table()
//   aim_fit_tree()                StemFinder::find()
//   aim_fit_slice()               StemFinder::fit_slice()
//   aim_reassign_ids()            InventoryMerger::reassign()
//   aim_resolve_multi_matching()  InventoryMerger::resolve()
//
// Differences with the R code are flagged with "R difference" in the comments.

namespace arbor::inventory {

using utils::fitting::Vec3;

namespace {

// false for "no tree" (0, -1 or R's NA)
inline bool is_tree(int id) { return id > 0; }

// R unique() keeping the order of first appearance
template <typename T>
void push_unique(std::vector<T>& v, T x)
{
  if (std::find(v.begin(), v.end(), x) == v.end()) v.push_back(x);
}

void require(const PointCloud& pc, const char* fun, bool hag, bool foliage)
{
  auto missing = [fun](const char* attr) { return std::runtime_error(std::string(fun) + ": point cloud is missing required '" + attr + "' attribute."); };
  if (!pc.has_treeid())           throw missing("treeID");
  if (hag && !pc.has_hag())       throw missing("hag");
  if (foliage && !pc.has_foliage()) throw missing("foliage");
}

void remap_ids(PointCloud& pc, const std::unordered_map<int, int>& map)
{
  if (map.empty()) return;
  for (size_t i = 0; i < pc.size(); ++i)
  {
    auto it = map.find(pc.get_treeid(i));
    if (it != map.end()) pc.set_treeid(i, it->second);
  }
}

// ---------------------------------------------------------------------------
// Stem search around an inventory tree
// ---------------------------------------------------------------------------

// 2D KD-tree on the points that can be useful (hag <= max(POM) + 2.5). The R code scans the
// whole point cloud for every inventory tree.
struct XYAdaptor
{
  const PointCloud& pc;
  std::vector<unsigned int> idx; // indices of the indexed points in pc

  inline size_t kdtree_get_point_count() const { return idx.size(); }
  inline double kdtree_get_pt(size_t i, size_t dim) const { return dim == 0 ? pc.get_x(idx[i]) : pc.get_y(idx[i]); }
  template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
};

using XYTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, XYAdaptor>, XYAdaptor, 2>;

// A point of the neighbourhood of an inventory tree, copied because it is accessed many times
struct NPoint
{
  double x, y, z, hag;
  int id;
};

class StemFinder
{
public:
  StemFinder(const PointCloud& pc, double max_hag, const MergerParameters& params) :
    pc_(pc), params_(params), adaptor_{pc, {}}, kdtree_(2, adaptor_, nanoflann::KDTreeSingleIndexAdaptorParams(10, nanoflann::KDTreeSingleIndexAdaptorFlags::SkipInitialBuildIndex))
  {
    for (size_t i = 0; i < pc.size(); ++i)
      if (pc.get_hag(i) <= max_hag) adaptor_.idx.push_back(static_cast<unsigned int>(i));
    kdtree_.buildIndex();
  }

  // The stem is searched at several heights around the POM in a disc of 1.5 radius. If nothing
  // is found, the disc grows to 3, 4.5 and 6 radii. The warning is returned instead of being
  // emitted because this runs in OpenMP threads, from which messages are discarded.
  StemMatch find(const InventoryTree& tree, std::string& warning) const
  {
    StemMatch res;
    res.target_id = tree.id;

    std::vector<NPoint> pts = neighbourhood(tree);

    std::vector<NPoint> candidates;
    double factor = params_.search_factor;
    for (int it = 0; it < params_.search_iterations && res.slices.empty(); ++it, factor += params_.search_factor)
    {
      double r2 = (tree.radius * factor) * (tree.radius * factor);
      candidates.clear();
      for (const NPoint& p : pts)
      {
        double dx = p.x - tree.x, dy = p.y - tree.y;
        if (dx*dx + dy*dy < r2) candidates.push_back(p);
      }

      for (double offset : params_.slice_offsets)
        if (auto slice = fit_slice(candidates, tree, offset)) res.slices.push_back(std::move(*slice));
    }

    // IDs to merge: union of the IDs of all slices. Seeds: inliers of all slices.
    if (!res.slices.empty())
    {
      res.status = StemMatch::Status::Fitted;
      for (const StemSlice& slice : res.slices)
      {
        for (int id : slice.tree_ids) push_unique(res.tree_ids, id);
        for (int k : slice.fit.inlier_indices) res.seeds.push_back(slice.points[k]);
      }
      return res;
    }

    // No stem found: if a single point cloud tree is present at the inventory position (+/- 0.5 m,
    // within 2 radii) it is used anyway. The seed is the inventory position itself.
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
      res.status = StemMatch::Status::SingleCandidate;
      res.tree_ids = ids;
      warning = "No tree detected for tree " + std::to_string(tree.id) + " but a single candidate ID was found and used.";
    }
    else
    {
      res.status = StemMatch::Status::NotDetected;
      warning = "No tree detected for tree " + std::to_string(tree.id);
    }

    return res;
  }

private:
  const PointCloud& pc_;
  const MergerParameters& params_;
  XYAdaptor adaptor_;
  XYTree kdtree_;

  // Points within 5 radii of the inventory position and below 3 m above the POM, clustered with
  // dbscan. Only the clusters that touch the ground are kept: this removes branches and crowns of
  // the neighbours hanging around the stem. Noise points are removed too.
  std::vector<NPoint> neighbourhood(const InventoryTree& tree) const
  {
    std::vector<nanoflann::ResultItem<XYTree::IndexType, double>> matches;
    double q[2] = {tree.x, tree.y};
    double r = tree.radius * params_.neighbourhood_factor;
    kdtree_.radiusSearch(q, r * r, matches, nanoflann::SearchParameters(0.0f, false));

    std::vector<unsigned int> idx;
    idx.reserve(matches.size());
    for (const auto& m : matches) idx.push_back(adaptor_.idx[m.first]);
    std::sort(idx.begin(), idx.end());

    std::vector<NPoint> pts;
    for (unsigned int i : idx)
    {
      double z = pc_.get_z(i);
      if (z < tree.z + params_.neighbourhood_above)
        pts.push_back({pc_.get_x(i), pc_.get_y(i), z, pc_.get_hag(i), pc_.get_treeid(i)});
    }

    if (pts.empty()) return pts;

    // Coordinates relative to the tree for numerical accuracy
    std::vector<point3> xyz;
    xyz.reserve(pts.size());
    for (const NPoint& p : pts) xyz.push_back({p.x - tree.x, p.y - tree.y, p.z - tree.z});

    // R difference: the vendored dbscan may attribute a border point reachable from two
    // clusters differently. This has no practical consequence.
    std::vector<std::vector<size_t>> clusters = dbscan(xyz, static_cast<float>(params_.cluster_eps), params_.cluster_min_pts);

    std::vector<NPoint> kept;
    for (const std::vector<size_t>& cluster : clusters)
    {
      double min_hag = std::numeric_limits<double>::infinity();
      for (size_t k : cluster) min_hag = std::min(min_hag, pts[k].hag);
      if (min_hag >= params_.cluster_max_base_hag) continue;
      for (size_t k : cluster) kept.push_back(pts[k]);
    }

    return kept;
  }

  // Fits the stem cross-section in a 30 cm thick horizontal slice centered at POM + offset
  std::optional<StemSlice> fit_slice(const std::vector<NPoint>& candidates, const InventoryTree& tree, double offset) const
  {
    // R difference (bug fix): the R code compares an absolute elevation with 1 m. The intent
    // is to skip slices lower than 1 m above ground.
    if (tree.pom + offset < params_.min_slice_height) return std::nullopt;

    double zc = tree.z + offset;
    double zmin = zc - params_.slice_half_thickness;
    double zmax = zc + params_.slice_half_thickness;

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

    // Not enough points, slice entirely close to the ground, or slice entirely much higher
    // than the POM (the stem does not reach the ground: not our tree)
    if (slice.size() < 3) return std::nullopt;
    if (max_hag < 0.5) return std::nullopt;
    if (min_hag > 2 * tree.pom) return std::nullopt;

    // Same as fit_circloid_cpp(xyz, complexity = 3): vertical axis through the centroid
    double mx = 0, my = 0;
    for (const NPoint* p : slice) { mx += p->x; my += p->y; }
    mx /= slice.size();
    my /= slice.size();

    utils::fitting::CrossSectionFitter fitter;
    fitter.set_axis({mx, my, 0}, {mx, my, 1});
    for (const NPoint* p : slice) fitter.add_point(p->x, p->y, p->z);

    utils::fitting::FittingResult res = fitter.fit(params_.fit_tolerance, utils::fitting::FitMode::Buttress);
    // A fit that covers less than half of the circle is not trusted to be a stem
    if (!res.success || res.arc_coverage_deg <= params_.min_arc_coverage) return std::nullopt;

    // R difference: points without tree (NA) are not included in the IDs to be merged
    StemSlice out;
    out.offset = offset;
    out.points.reserve(slice.size());
    for (const NPoint* p : slice) out.points.emplace_back(p->x, p->y, p->z);
    for (int k : res.inlier_indices)
      if (is_tree(slice[k]->id)) push_unique(out.tree_ids, slice[k]->id);
    out.fit = std::move(res);
    return out;
  }
};

// The matching table gives rules treeID -> targetID. The R code works in 3 phases on the IDs:
//   Phase 0: pre-existing IDs that are equal to a target (and that are not themselves a
//            source) are moved to fresh IDs, otherwise they would merge with the points
//            relabelled to this target.
//   Phase 1: source -> temporary ID
//   Phase 2: temporary ID -> target
// The final ID of a point only depends on its original ID, so the 3 phases are collapsed
// into a single lookup table old ID -> new ID. Fresh IDs are the same as in R.
void reassign_ids(PointCloud& pc, const MatchingTable& table, bool warn)
{
  require(pc, "reassign", false, false);

  for (const Match& m : table)
    if (m.target_id <= 0) throw std::invalid_argument("reassign: invalid target ID " + std::to_string(m.target_id) + " in the matching table.");

  std::map<int, std::vector<int>> conflicts = find_conflicts(table);
  if (warn && !conflicts.empty())
    logging::warning("Matching table contains duplicated 'treeID': use resolve().");

  // Rules source -> target. Rows without tree, conflicting rows (handled by resolve()) and
  // identity rows are removed. Identity rules must not be considered as handled IDs, otherwise
  // they are not protected from the points reassigned to them.
  std::vector<int> src, tgt;
  std::unordered_set<int> seen;
  for (const Match& m : table)
  {
    if (!is_tree(m.tree_id) || m.tree_id == m.target_id) continue;
    if (conflicts.count(m.tree_id) || !seen.insert(m.tree_id).second) continue;
    src.push_back(m.tree_id);
    tgt.push_back(m.target_id);
  }

  if (src.empty()) return;

  // Targets that are not sources collide with the pre-existing IDs of the point cloud: these
  // IDs are moved to new IDs. A target that is also a source is freed when its source moves.
  std::unordered_set<int> sources(src.begin(), src.end());
  std::unordered_set<int> colliding;
  for (int t : tgt)
    if (!sources.count(t)) colliding.insert(t);

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

  // Fresh IDs: max_val + n_src + k (the first n_src IDs are the temporary IDs of phase 1 in R)
  std::unordered_map<int, int> map;
  for (size_t k = 0; k < src.size(); ++k) map[src[k]] = tgt[k];
  for (size_t k = 0; k < bump.size(); ++k) map[bump[k]] = max_val + static_cast<int>(src.size() + k + 1);

  remap_ids(pc, map);
}

} // namespace

// ---------------------------------------------------------------------------
// Matching table
// ---------------------------------------------------------------------------

MatchingTable matching_table(const std::vector<StemMatch>& matches)
{
  MatchingTable table;
  for (const StemMatch& m : matches)
  {
    if (m.tree_ids.empty())
    {
      table.push_back({0, m.target_id, 0});
      continue;
    }

    for (int id : m.tree_ids)
      table.push_back({id, m.target_id, m.nfit()});
  }
  return table;
}

// Typically one arbor tree that is in reality two or more inventory trees.
// R difference: R counts duplicated rows, here identical (treeID, targetID) rows given
// twice are not a conflict.
std::map<int, std::vector<int>> find_conflicts(const MatchingTable& table)
{
  std::map<int, std::vector<int>> targets;
  for (const Match& m : table)
    if (is_tree(m.tree_id)) push_unique(targets[m.tree_id], m.target_id);

  std::map<int, std::vector<int>> conflicts;
  for (auto& [id, t] : targets)
    if (t.size() > 1) conflicts.emplace(id, std::move(t));
  return conflicts;
}

// ---------------------------------------------------------------------------
// Report
// ---------------------------------------------------------------------------

std::vector<int> MergeReport::undetected() const
{
  std::vector<int> ids;
  for (const StemMatch& m : matches)
    if (m.tree_ids.empty()) ids.push_back(m.target_id);
  return ids;
}

std::string MergeReport::to_string() const
{
  std::vector<int> missing = undetected();
  std::map<int, std::vector<int>> dup = conflicts();

  std::ostringstream os;
  os << "Inventory merging: " << matches.size() << " inventory trees | "
     << matches.size() - missing.size() << " matched | "
     << missing.size() << " undetected | "
     << dup.size() << " conflicts | "
     << disentangled.size() << " IDs renumbered";

  if (!missing.empty())
  {
    os << "\nUndetected:";
    for (int id : missing) os << " " << id;
  }

  if (!dup.empty())
  {
    os << "\nConflicts (re-segmented):";
    for (const auto& [id, t] : dup)
    {
      os << " " << id << "->";
      for (size_t k = 0; k < t.size(); ++k) os << (k ? "," : "") << t[k];
    }
  }

  return os.str();
}

// ---------------------------------------------------------------------------
// InventoryMerger
// ---------------------------------------------------------------------------

InventoryMerger::InventoryMerger(const Inventory& inventory, MergerParameters params) : inventory_(inventory), params_(std::move(params))
{
}

// If a point cloud tree has an ID that also exists in the inventory (e.g. arbor's tree 54 and
// inventory's tree 54, almost surely two different trees), the point cloud tree gets a new ID.
// After that, point cloud IDs and inventory IDs never overlap and reassign() cannot confuse them.
std::unordered_map<int, int> InventoryMerger::disentangle(PointCloud& pc) const
{
  require(pc, "disentangle", false, false);

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

  // R difference (bug fix): new IDs start after the max of both the point cloud and the
  // inventory IDs. In R a new ID could be equal to an inventory ID.
  for (const InventoryTree& t : inventory_) max_id = std::max(max_id, t.id);

  std::unordered_map<int, int> map;
  int next = max_id + 1;
  for (int id : ids)
    if (inventory_.contains(id)) map[id] = next++;

  if (map.empty())
    logging::info("No matching treeIDs found between the point cloud and the inventory.");

  remap_ids(pc, map);
  return map;
}

std::vector<StemMatch> InventoryMerger::match(const PointCloud& pc) const
{
  require(pc, "match", true, false);

  for (size_t i = 0; i < pc.size(); ++i)
  {
    if (inventory_.contains(pc.get_treeid(i)))
      throw std::runtime_error("Some tree IDs in the point cloud are also in the inventory tree IDs. This can cause matching issues. Call disentangle() first.");
  }

  std::vector<const InventoryTree*> trees;
  trees.reserve(inventory_.size());
  for (const InventoryTree& t : inventory_) trees.push_back(&t);
  return fit(pc, trees);
}

std::vector<StemMatch> InventoryMerger::match(const PointCloud& pc, const std::vector<int>& ids) const
{
  require(pc, "match", true, false);

  std::vector<const InventoryTree*> trees;
  trees.reserve(ids.size());
  for (int id : ids) trees.push_back(&inventory_.at(id));
  return fit(pc, trees);
}

// The trees are independent and are processed in parallel. Warnings are collected per tree and
// emitted after the loop, in inventory order. An exception in a thread stops the loop and is
// rethrown afterwards.
std::vector<StemMatch> InventoryMerger::fit(const PointCloud& pc, const std::vector<const InventoryTree*>& trees) const
{
  std::vector<StemMatch> matches(trees.size());
  if (trees.empty()) return matches;

  // Points above max(POM) + 2.5 are never useful
  double max_pom = 0;
  for (const InventoryTree& t : inventory_) max_pom = std::max(max_pom, t.pom);
  StemFinder finder(pc, max_pom + params_.max_hag_above_pom, params_);

  std::vector<std::string> warnings(trees.size());
  std::exception_ptr error = nullptr;
  std::atomic<bool> failed{false};

  Counter counter("Fitting inventory trees", trees.size());

  #pragma omp parallel for schedule(dynamic)
  for (size_t k = 0; k < trees.size(); ++k)
  {
    if (counter.cancelled() || failed.load(std::memory_order_relaxed)) continue;

    try
    {
      matches[k] = finder.find(*trees[k], warnings[k]);
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

  return matches;
}

void InventoryMerger::reassign(PointCloud& pc, const MatchingTable& table)
{
  reassign_ids(pc, table, true);
}

// When one arbor tree (e.g. 123) is matched to several inventory trees (e.g. 12 and 13), the
// points of tree 123 are re-segmented with segment_instance() using the stem points found for
// 12 and 13 as seeds. Each point then goes to the closest seed in the graph sense, which splits
// the tree. All the conflicting trees are re-segmented together in a single call.
void InventoryMerger::resolve(PointCloud& pc, const MatchingTable& table, const std::vector<StemMatch>& matches) const
{
  std::map<int, std::vector<int>> conflicts = find_conflicts(table);
  if (conflicts.empty()) return;

  require(pc, "resolve", true, true);

  std::unordered_map<int, const StemMatch*> match_of;
  for (const StemMatch& m : matches) match_of[m.target_id] = &m;

  // Inventory trees involved without a match are fitted again.
  // R difference: the fits use the same point cloud filter as the matching step (hag <= max(POM)
  // + 2.5) instead of the full point cloud, which is faster and consistent with the matching.
  std::vector<const InventoryTree*> missing;
  for (const auto& [id, targets] : conflicts)
  {
    for (int target : targets)
    {
      if (match_of.count(target)) continue;
      match_of[target] = nullptr;
      missing.push_back(&inventory_.at(target));
    }
  }

  std::vector<StemMatch> refitted = fit(pc, missing);
  for (const StemMatch& m : refitted) match_of[m.target_id] = &m;

  // Seeds: stem inliers of each inventory tree involved
  std::vector<std::pair<Vec3, int>> seed_pts;
  for (const auto& [id, targets] : conflicts)
    for (int target : targets)
      for (const Vec3& p : match_of[target]->seeds) seed_pts.emplace_back(p, target);

  // Points of the conflicting trees. 'points' is the index of each point of the subset in pc.
  std::vector<unsigned int> points;
  for (size_t i = 0; i < pc.size(); ++i)
    if (conflicts.count(pc.get_treeid(i))) points.push_back(static_cast<unsigned int>(i));

  if (points.empty() || seed_pts.empty()) return;

  // segment_instance() reorders the points (partition above/below cut_above_ground). The
  // original index is stored in the 'passage' attribute, which segment_instance() moves with
  // the points and does not use.
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

  // R difference: the oversegmentation solver is disabled, it would merge the trees we are splitting
  settings::ArborParameters params;
  params.instance.oversegmentation_solver_enabled = false;
  segment::segment_instance(sub, seeds, params);

  // R difference: points not reached by the segmentation (NA) or below cut_above_ground (0)
  // keep their original ID instead of losing it. true_size() includes the points moved
  // beyond size() by the partition.
  for (size_t k = 0; k < sub.true_size(); ++k)
  {
    int id = sub.get_treeid(k);
    if (is_tree(id)) pc.set_treeid(points[sub.get_passage(k)], id);
  }
}

// The fits of the matching step are reused to resolve the conflicts instead of being recomputed
MergeReport InventoryMerger::run(PointCloud& pc) const
{
  require(pc, "merge", true, false);
  if (inventory_.empty()) throw std::runtime_error("merge: inventory is empty.");

  MergeReport report;

  Task task("Inventory merging", {1, 20, 1, 10});

  task.step("Disentangling IDs");
  report.disentangled = disentangle(pc);

  task.step("Matching inventory trees");
  report.matches = match(pc);
  MatchingTable table = report.table();

  task.step("Reassigning IDs");
  reassign_ids(pc, table, false); // conflicts are expected and handled below

  task.step("Resolving multi-matching");
  resolve(pc, table, report.matches);

  logging::info(report.to_string());

  return report;
}

} // namespace arbor::inventory

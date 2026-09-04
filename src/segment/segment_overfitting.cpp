/**
 * @file segment_overfitting.cpp
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

#include "segment_overfitting.h"
#include "dbscan/dbscan.hpp"
#include "fitting_orbicular.h"
#include "myomp.h"
#include "services.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace arbor::segment {

namespace {

inline int layer_count(const MergerConfig& cfg)
{
  return static_cast<int>(std::llround((cfg.layer_end - cfg.layer_start) / cfg.layer_step)) + 1;
}

inline int layer_of(double hag, const MergerConfig& cfg, int nlayers)
{
  if (hag < cfg.layer_start) return -1;
  const int i = static_cast<int>(std::floor((hag - cfg.layer_start) / cfg.layer_step + 1e-9));
  if (i < 0 || i >= nlayers) return -1;
  return i;
}

// Packs a sorted (ref, move) treeID pair into a single 64-bit key so pairs
// can be deduplicated/accumulated in a hash map. Callers must ensure
// ref < move (true here since `tid` is sorted before matches are built).
inline int64_t make_match_key(int ref, int move)
{
  return (static_cast<int64_t>(ref) << 32) | static_cast<uint32_t>(move);
}

// Disjoint-set helper used to collapse chains of matches (e.g. 1-2, 2-3)
// onto a single canonical treeID before remapping the point cloud, so a
// merge never targets an id that a previous merge already remapped away.
class UnionFind
{
public:
  int find(int id)
  {
    auto it = parent_.find(id);
    if (it == parent_.end() || it->second == id) return id;
    int root = find(it->second);
    it->second = root; // path compression
    return root;
  }

  // Unions the two ids, keeping the smaller value as the canonical root.
  void unite(int a, int b)
  {
    int ra = find(a);
    int rb = find(b);
    if (ra == rb) return;
    int keep    = std::min(ra, rb);
    int discard = std::max(ra, rb);
    parent_[discard] = keep;
  }

  const std::unordered_map<int, int>& parent() const { return parent_; }

private:
  std::unordered_map<int, int> parent_;
};

} // namespace

OverSegmentationResolver::OverSegmentationResolver(MergerConfig config) : config_(config)
{
}

InstanceMatchResult OverSegmentationResolver::detect(const PointCloud& cloud) const
{
  if (!cloud.has_hag())
    throw std::runtime_error("OverSegmentationResolver::detect: point cloud has no HAG attribute");
  if (!cloud.has_treeid())
    throw std::runtime_error("OverSegmentationResolver::detect: point cloud has no treeID attribute");

  const size_t npoints = cloud.size();
  const int nlayers = layer_count(config_);

  ServiceLocator::logger()("Overfitting enhancement module");

  // Parallel point bucketing into HAG layers. Each thread accumulates into
  // its own local buckets and merges them once at the end, instead of
  // locking per point, so the critical section only ever runs `nlayers`
  // times per thread rather than `npoints` times.
  std::vector<std::vector<size_t>> layer_indices(nlayers);
  for (auto& bucket : layer_indices)
    bucket.reserve(npoints / static_cast<size_t>(nlayers)); // rough heuristic, avoids most reallocations

  #pragma omp parallel
  {
    std::vector<std::vector<size_t>> local_buckets(nlayers);

    #pragma omp for
    for (size_t point_idx = 0; point_idx < npoints; ++point_idx)
    {
      const int layer_idx = layer_of(cloud.get_hag(point_idx), config_, nlayers);
      if (layer_idx >= 0)
      {
        local_buckets[layer_idx].push_back(point_idx);
      }
    }

    #pragma omp critical
    {
      for (int layer_idx = 0; layer_idx < nlayers; ++layer_idx)
      {
        layer_indices[layer_idx].insert(layer_indices[layer_idx].end(), local_buckets[layer_idx].begin(), local_buckets[layer_idx].end());
      }
    }
  }

  // Parallel layer processing
  std::vector<InstanceMatchResult> layer_results(nlayers);

  auto pb = ServiceLocator::make_progress(nlayers, "Searching overfitting");
  std::atomic<bool> abort(false);

  // For each layer
  #pragma omp parallel for schedule(dynamic)
  for (int layer_idx = 0; layer_idx < nlayers; ++layer_idx)
  {
    if (abort.load(std::memory_order_relaxed)) continue;
    if(pb->check_interrupt()) abort = true;
    pb->tick();

    const std::vector<size_t>& layer_point_indices = layer_indices[layer_idx];

    if (layer_point_indices.size() < static_cast<size_t>(config_.dbscan_min_pts))
    {
      continue;
    }

    // DBSCAN on this layer
    std::vector<point3> pts;
    pts.reserve(layer_point_indices.size());
    for (size_t point_idx : layer_point_indices)
    {
      pts.push_back(point3{cloud.get_x(point_idx), cloud.get_y(point_idx), cloud.get_z(point_idx)});
    }

    std::vector<std::vector<size_t>> clusters = dbscan(pts, config_.dbscan_eps, 1);

    // For each clusters
    for (const auto& cluster : clusters)
    {
      if (cluster.size() < static_cast<size_t>(config_.dbscan_min_pts))
      {
        continue;
      }

      // Centroid calculation to center on 0,0 and improve floating point accuracy stability
      std::vector<size_t> global_indices;
      global_indices.reserve(cluster.size());
      double mean_x = 0.0;
      double mean_y = 0.0;
      double mean_hag = 0.0; // Track cluster HAG

      for (size_t local_idx : cluster)
      {
        const size_t point_idx = layer_point_indices[local_idx];
        global_indices.push_back(point_idx);
        mean_x += cloud.get_x(point_idx);
        mean_y += cloud.get_y(point_idx);
        mean_hag += cloud.get_hag(point_idx);
      }

      const double inv_size = 1.0 / static_cast<double>(cluster.size());
      mean_x *= inv_size;
      mean_y *= inv_size;
      mean_hag *= inv_size;

      // Fit cross-section on the cluster
      utils::fitting::CrossSectionFitter fitter;
      fitter.set_axis(utils::fitting::Vec3(mean_x, mean_y, 0.0), utils::fitting::Vec3(mean_x, mean_y, 1.0));
      for (size_t point_idx : global_indices)
      {
        fitter.add_point(cloud.get_x(point_idx), cloud.get_y(point_idx), cloud.get_z(point_idx));
      }

      utils::fitting::FittingResult fit = fitter.fit(config_.fit_tolerance, utils::fitting::FitMode::Buttress);

      // Rejection filters. Reject circle that are not of good enough quality
      if (!fit.success ||
          fit.arc_coverage_deg < config_.min_arc_degree ||
          fit.inlier_percentage < config_.min_inlier_pct ||
          fit.interior_percentage > config_.max_interior_pct ||
          fit.radius > config_.max_radius ||
          fit.radius < config_.min_radius)
      {
        continue;
      }

      // For tree with < 10 cm radius. No need to search too high. This is prone to overfitting. Problem
      // Arise mainly for big trees.
      if (fit.radius < 0.2  && mean_hag > 6.0) continue;
      if (fit.radius < 0.1  && mean_hag > 4.0) continue;
      if (fit.radius < 0.05 && mean_hag > 2.0) continue;

      // Count inlier points per treeID supporting the circle. If id_count = 1 then the circle has been
      // fit on a single tree. The instance segmentation was good. If id_count > 1 then a circle has been
      // fit on points with more than 1 ID. This might be an over-segmentation.
      //
      // A plain vector with linear lookup is used instead of a hash map: a
      // circle almost never supports more than a handful of distinct
      // treeIDs (anything with only 1 is discarded below), so the hashing
      // and allocation overhead of an unordered_map isn't worth paying for
      // every single one of the (potentially very many) accepted circles.
      std::vector<std::pair<int, size_t>> id_counts;
      id_counts.reserve(4);
      for (int local_inlier : fit.inlier_indices)
      {
        if (local_inlier < 0 || static_cast<size_t>(local_inlier) >= cluster.size())
          continue;

        const int tree_id = cloud.get_treeid(global_indices[local_inlier]);
        auto it = std::find_if(id_counts.begin(), id_counts.end(),
                                [tree_id](const auto& kv) { return kv.first == tree_id; });
        if (it == id_counts.end())
          id_counts.emplace_back(tree_id, 1);
        else
          ++it->second;
      }

      // id_count = 1 then a circle has been fit on points with more than 1 ID. No problem here
      if (id_counts.size() < 2)
      {
        continue;
      }

      std::vector<int> tree_ids;
      tree_ids.reserve(id_counts.size());
      for (const auto& kv : id_counts) tree_ids.push_back(kv.first);
      std::sort(tree_ids.begin(), tree_ids.end());

      auto points_for = [&id_counts](int tree_id) -> size_t
      {
        auto it = std::find_if(id_counts.begin(), id_counts.end(),
                                [tree_id](const auto& kv) { return kv.first == tree_id; });
        return it != id_counts.end() ? it->second : 0;
      };

      // Record circle and match pairs into thread-isolated storage
      layer_results[layer_idx].circles.push_back(CircleRecord{fit.center.x, fit.center.y, fit.center.z, fit.radius, tree_ids.size()});

      for (size_t k = 1; k < tree_ids.size(); ++k)
      {
        layer_results[layer_idx].matches.push_back(MatchRecord{tree_ids[0], tree_ids[k], 1, points_for(tree_ids[0]), points_for(tree_ids[k])});
      }
    }
  }

  if (abort.load()) throw std::runtime_error("Computation aborted");

  pb->finalize();

  // Merge of per-layer results into final structure
  InstanceMatchResult result;
  size_t total_circles = 0;
  size_t total_matches = 0;

  for (const auto& res : layer_results)
  {
    total_circles += res.circles.size();
    total_matches += res.matches.size();
  }

  result.circles.reserve(total_circles);

  // Collapse duplicate (ref, move) pairs into a single row with an occurrence
  // count, instead of emitting one matching_table row per supporting circle.
  // tree_ids is sorted before matches are built, so ref < move always holds
  // for a given record; pack both into one 64-bit key for the hash map.
  //
  // MatchRecord is reused directly as the accumulator (rather than a
  // separate accumulator struct) since it already has exactly the fields
  // (n, ref_points, move_points) accumulation needs, plus ref/move to store
  // once up front -- this also avoids having to decode ref/move back out of
  // the key afterwards.
  std::unordered_map<int64_t, MatchRecord> match_counts;
  match_counts.reserve(total_matches);

  for (auto& res : layer_results)
  {
    result.circles.insert(
      result.circles.end(),
      std::make_move_iterator(res.circles.begin()),
      std::make_move_iterator(res.circles.end())
    );

    for (const MatchRecord& m : res.matches)
    {
      const int64_t key = make_match_key(m.ref, m.move);
      auto [it, inserted] = match_counts.try_emplace(key, MatchRecord{m.ref, m.move, 0, 0, 0});
      MatchRecord& acc = it->second;
      acc.n += m.n;
      acc.ref_points += m.ref_points;
      acc.move_points += m.move_points;
    }
  }

  result.matches.reserve(match_counts.size());
  for (auto& [key, acc] : match_counts)
  {
    result.matches.push_back(std::move(acc));
  }

  ServiceLocator::logger()(" Detection completed. Candidate found: " + std::to_string(result.matches.size()));

  return result;
}

void OverSegmentationResolver::merge(PointCloud& cloud, const std::vector<MatchRecord>& matches) const
{
  UnionFind uf;

  for (const MatchRecord& m : matches)
  {
    if (m.n <= config_.min_support) continue;
    if (config_.min_weight_ratio > 0.0 && m.weight_ratio() < config_.min_weight_ratio) continue;
    uf.unite(m.ref, m.move);
  }

  size_t merged_count = 0;

  for (const auto& entry : uf.parent())
  {
    const int id   = entry.first;
    const int root = uf.find(id);
    if (root != id)
    {
      cloud.remap_treeid(id, root);
      ServiceLocator::logger()(" Merged tree " + std::to_string(id) + " -> " + std::to_string(root));
      ++merged_count;
    }
  }

  ServiceLocator::logger()(" Merge completed. Total remapped trees: " + std::to_string(merged_count));
}

InstanceMatchResult OverSegmentationResolver::run(PointCloud& cloud) const
{
  InstanceMatchResult result = detect(cloud);
  merge(cloud, result.matches);
  return result;
}

} // namespace arbor::segment

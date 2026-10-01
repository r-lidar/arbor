/**
 * @file fitting_multicircle.cpp
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

#include "fitting_multicircle.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>

namespace arbor::utils::fitting {

MultiCircleFitter::MultiCircleFitter(int n_circles, int max_iterations, double early_exit_ratio,
                                     unsigned seed, double max_overlap_ratio, double min_inlier_gain_ratio)
  : m_n_circles(n_circles),
    m_max_iterations(max_iterations),
    m_early_exit_ratio(early_exit_ratio),
    m_max_overlap_ratio(max_overlap_ratio),
    m_min_inlier_gain_ratio(min_inlier_gain_ratio),
    m_rng(seed)
{
}

FittingResult MultiCircleFitter::fit(const std::vector<Vec3>& points, double tolerance)
{
  FittingResult result;
  result.shape_type = "multi_circle_" + std::to_string(m_n_circles);

  if (static_cast<int>(points.size()) < 3 * m_n_circles)
  {
    result.success = false;
    return result;
  }

  double zsum = 0.0;
  for (const auto& v : points) zsum += v.z;
  m_zmean = zsum / static_cast<double>(points.size());

  // --- Sequential RANSAC ---
  // Fit the best circle on the remaining unclaimed points, then remove those
  // inliers before fitting the next circle. The first circle fitted here is,
  // by construction, the best single circle for the *whole* point set, it
  // becomes the baseline for the inlier-gain guard below.
  std::vector<int> remaining(points.size());
  std::iota(remaining.begin(), remaining.end(), 0);

  std::vector<CircleModel>       fitted_circles;
  std::vector<std::vector<int>>  per_circle_inliers;

  for (int c = 0; c < m_n_circles; ++c)
  {
    if (static_cast<int>(remaining.size()) < 3) break;

    CircleModel circle = fit_circle_ransac(points, remaining, tolerance);
    if (!circle.valid) break;

    std::vector<int> inliers = find_inliers(points, remaining, circle, tolerance);
    if (inliers.empty()) break;

    fitted_circles.push_back(circle);
    per_circle_inliers.push_back(inliers);

    std::set<int> inlier_set(inliers.begin(), inliers.end());
    std::vector<int> new_remaining;
    new_remaining.reserve(remaining.size() - inliers.size());
    for (int idx : remaining)
    {
        if (!inlier_set.count(idx))
            new_remaining.push_back(idx);
    }
    remaining = std::move(new_remaining);
  }

  if (static_cast<int>(fitted_circles.size()) < m_n_circles)
  {
    result.success = false;
    return result;
  }

  // --- Overlap guard ---
  // Reject if any two fitted circles overlap too much. Overlap is measured
  // as the intersection area relative to the area of the smaller of the two
  // circles, so a small circle nested inside a big one is caught even when
  // their centers are far apart relative to the big circle's radius.
  for (int i = 0; i < static_cast<int>(fitted_circles.size()); ++i)
  {
    for (int j = i + 1; j < static_cast<int>(fitted_circles.size()); ++j)
    {
      if (circle_overlap_ratio(fitted_circles[i], fitted_circles[j]) > m_max_overlap_ratio)
      {
        result.success = false;
        return result;
      }
    }
  }

  // --- Inlier-gain guard ---
  // The extra circle(s) must explain meaningfully more of the point cloud
  // than the single best circle already does. Without this check, a second
  // circle fitted on the leftover outliers from the first one is almost
  // always "found" even for genuinely single-stem data, since RANSAC will
  // happily fit a circle to a handful of scattered outlier points. Require
  // the combined inlier count across all circles to exceed the first
  // (best single-circle) fit's inlier count by at least the configured
  // gain ratio (e.g. 1.3 = +30%).
  if (fitted_circles.size() > 1)
  {
    const double baseline_inliers = static_cast<double>(per_circle_inliers.front().size());
    double total_inliers = 0.0;
    for (const auto& inl : per_circle_inliers) total_inliers += static_cast<double>(inl.size());

    if (total_inliers < m_min_inlier_gain_ratio * baseline_inliers)
    {
      result.success = false;
      return result;
    }
  }

  // --- Aggregate metrics ---

  // Arc coverage: max across all circles (less conservative).
  double min_arc = std::numeric_limits<double>::max();
  for (int i = 0; i < static_cast<int>(fitted_circles.size()); ++i)
  {
    const Vec3 c3d = {fitted_circles[i].cx, fitted_circles[i].cy, m_zmean};
    min_arc = std::max(min_arc, calculate_arc_coverage(points, per_circle_inliers[i], c3d));
  }

  // Interior percentage: a point counts if it falls strictly inside any circle.
  int interiors = 0;
  for (const auto& p : points)
  {
    for (const auto& c : fitted_circles)
    {
      const double dx = p.x - c.cx;
      const double dy = p.y - c.cy;
      if (std::sqrt(dx*dx + dy*dy) < c.radius - tolerance) { ++interiors; break; }
    }
  }

  // Combined center: area-weighted centroid.
  double total_area = 0.0, cx_w = 0.0, cy_w = 0.0;
  for (const auto& c : fitted_circles)
  {
    const double area = c.radius * c.radius;
    total_area += area;
    cx_w += c.cx * area;
    cy_w += c.cy * area;
  }
  cx_w /= total_area;
  cy_w /= total_area;

  // Equivalent radius: preserves total cross-section area.
  double r_equiv = 0.0;
  for (const auto& c : fitted_circles) r_equiv += c.radius * c.radius;
  r_equiv = std::sqrt(r_equiv);

  // Collect all inliers
  std::vector<int> all_inliers;
  for (const auto& inl : per_circle_inliers)
    all_inliers.insert(all_inliers.end(), inl.begin(), inl.end());

  // --- Fill result ---
  result.success             = true;
  result.center              = {cx_w, cy_w, m_zmean};
  result.radius              = static_cast<float>(r_equiv);
  result.inlier_indices      = all_inliers;
  result.inlier_percentage   = 100.0f * static_cast<float>(all_inliers.size()) / static_cast<float>(points.size());
  result.interior_percentage = 100.0f * static_cast<float>(interiors)          / static_cast<float>(points.size());
  result.arc_coverage_deg    = static_cast<float>(min_arc);

  // Parameters packed as [cx₁, cy₁, r₁,  cx₂, cy₂, r₂, …]
  for (const auto& c : fitted_circles)
  {
    result.parameters.push_back(c.cx);
    result.parameters.push_back(c.cy);
    result.parameters.push_back(c.radius);
  }

  // Contour: full outline for each circle, concatenated
  for (const auto& c : fitted_circles)
  {
    for (int deg = 0; deg <= 360; deg += 2)
    {
      const double t = deg * M_PI / 180.0;
      result.contour.push_back({c.cx + c.radius * std::cos(t),
                               c.cy + c.radius * std::sin(t),
                               m_zmean});
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

CircleModel MultiCircleFitter::fit_circle_ransac(
    const std::vector<Vec3>& points,
    const std::vector<int>&  indices,
    double                   tolerance) const
{
  CircleModel best;
  const int n = static_cast<int>(indices.size());
  if (n < 3) return best;

  //std::uniform_int_distribution<int> dist(0, n - 1);
  deterministic_uniform_int_distribution dist(0, n - 1);
  int max_inliers           = 0;
  const int early_exit_thr  = static_cast<int>(m_early_exit_ratio * n);

  for (int iter = 0; iter < m_max_iterations; ++iter)
  {
    int li1 = dist(m_rng);
    int li2; do { li2 = dist(m_rng); } while (li2 == li1);
    int li3; do { li3 = dist(m_rng); } while (li3 == li1 || li3 == li2);

    const Vec3& p1 = points[indices[li1]];
    const Vec3& p2 = points[indices[li2]];
    const Vec3& p3 = points[indices[li3]];

    CircleModel circle = fit_circle_on_3_points(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
    if (!circle.valid || circle.radius <= 0.0 || std::isnan(circle.radius)) continue;

    int inlier_count = 0;
    for (int idx : indices)
    {
      if (point_to_circle_distance(points[idx].x, points[idx].y, circle) < tolerance)
        ++inlier_count;
    }

    if (inlier_count > max_inliers)
    {
        max_inliers = inlier_count;
        best = circle;
        if (max_inliers >= early_exit_thr) break;
    }
  }

  if (max_inliers > 0) best.valid = true;
  return best;
}

std::vector<int> MultiCircleFitter::find_inliers(
    const std::vector<Vec3>& points,
    const std::vector<int>&  candidates,
    const CircleModel&       circle,
    double                   tolerance) const
{
  std::vector<int> inliers;
  for (int idx : candidates)
  {
    if (point_to_circle_distance(points[idx].x, points[idx].y, circle) <= tolerance)
      inliers.push_back(idx);
  }
  return inliers;
}

} // namespace arbor::utils::fitting

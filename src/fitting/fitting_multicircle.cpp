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

#include "fitting.h"
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include <set>

namespace arbor::utils::fitting {

FittingMultiCircle::FittingMultiCircle(int n_circles, int max_iterations, double early_exit_ratio, unsigned seed)
  : m_n_circles(n_circles),
    m_max_iterations(max_iterations),
    m_early_exit_ratio(early_exit_ratio),
    m_rng(seed)
{
}

FittingResult FittingMultiCircle::fit(const std::vector<Vec3>& points, double tolerance)
{
  FittingResult result;
  result.shape_type = "multi_circle_" + std::to_string(m_n_circles);

  // Need at least 3 points per circle to attempt a fit
  if (static_cast<int>(points.size()) < 3 * m_n_circles)
  {
    result.success = false;
    return result;
  }

  double zsum = 0.0;
  for (const auto& v : points) zsum += v.z;
  m_zmean = zsum / static_cast<double>(points.size());

  // --- Sequential RANSAC ---
  // Each iteration fits the best circle on the remaining unclaimed points,
  // then removes those inliers before the next circle is attempted.
  std::vector<int> remaining(points.size());
  std::iota(remaining.begin(), remaining.end(), 0);

  std::vector<CircleParams>       fitted_circles;
  std::vector<std::vector<int>>   per_circle_inliers;

  for (int c = 0; c < m_n_circles; ++c)
  {
    if (static_cast<int>(remaining.size()) < 3) break;

    CircleParams circle = fit_circle_ransac(points, remaining, tolerance);
    if (!circle.valid) break;

    std::vector<int> inliers = find_inliers(points, remaining, circle, tolerance);
    if (inliers.empty()) break;

    fitted_circles.push_back(circle);
    per_circle_inliers.push_back(inliers);

    // Remove the claimed inliers from the candidate pool
    std::set<int> inlier_set(inliers.begin(), inliers.end());
    std::vector<int> new_remaining;
    new_remaining.reserve(remaining.size() - inliers.size());
    for (int idx : remaining)
      if (!inlier_set.count(idx))
        new_remaining.push_back(idx);
      remaining = std::move(new_remaining);
  }

  // Only succeed if we actually found all requested circles
  if (static_cast<int>(fitted_circles.size()) < m_n_circles)
  {
    result.success = false;
    return result;
  }

  // --- Separation guard ---
  // Reject if any two circles are so close their centers overlap one another,
  // which indicates they are redundantly fitting the same cluster.
  for (int i = 0; i < static_cast<int>(fitted_circles.size()); ++i)
  {
    for (int j = i + 1; j < static_cast<int>(fitted_circles.size()); ++j)
    {
      const auto& ci = fitted_circles[i];
      const auto& cj = fitted_circles[j];
      double dx   = ci.cx - cj.cx;
      double dy   = ci.cy - cj.cy;
      double dist = std::sqrt(dx * dx + dy * dy);
      // Circles are considered the same if their centers are within 30 % of the
      // combined radii — far tighter than two adjacent stems would ever be.
      if (dist < 0.3 * (ci.radius + cj.radius))
      {
        result.success = false;
        return result;
      }
    }
  }

  // --- Aggregate metrics ---

  // Arc coverage: use the minimum across all circles (most conservative).
  // If any individual circle has poor coverage, the full fit is unreliable.
  double min_arc = std::numeric_limits<double>::max();
  for (int i = 0; i < static_cast<int>(fitted_circles.size()); ++i)
  {
    Vec3 c3d = { fitted_circles[i].cx, fitted_circles[i].cy, m_zmean };
    double arc = calculate_arc_coverage(points, per_circle_inliers[i], c3d);
    min_arc = std::min(min_arc, arc);
  }

  // Insider percentage: a point counts if it falls strictly inside *any* circle.
  int insiders = 0;
  for (const auto& p : points)
  {
    for (const auto& c : fitted_circles)
    {
      double dx = p.x - c.cx;
      double dy = p.y - c.cy;
      if (std::sqrt(dx * dx + dy * dy) < c.radius - tolerance)
      {
        ++insiders;
        break; // each point counted once regardless of how many circles it's inside
      }
    }
  }

  // Combined center: area-weighted centroid (radius² ∝ cross-section area).
  double total_area = 0.0;
  double cx_w = 0.0, cy_w = 0.0;
  for (const auto& c : fitted_circles)
  {
    double area = c.radius * c.radius;
    total_area += area;
    cx_w += c.cx * area;
    cy_w += c.cy * area;
  }
  cx_w /= total_area;
  cy_w /= total_area;

  // Equivalent radius: sqrt(Σ rᵢ²) — preserves total cross-section area.
  double r_equiv = 0.0;
  for (const auto& c : fitted_circles)
    r_equiv += c.radius * c.radius;
  r_equiv = std::sqrt(r_equiv);

  // Collect all inliers
  std::vector<int> all_inliers;
  for (const auto& inl : per_circle_inliers)
    all_inliers.insert(all_inliers.end(), inl.begin(), inl.end());

  // --- Fill result ---
  result.success            = true;
  result.center             = { cx_w, cy_w, m_zmean };
  result.radius             = static_cast<float>(r_equiv);
  result.inlier_indices     = all_inliers;
  result.inlier_percentage  = 100.0 * static_cast<double>(all_inliers.size()) / points.size();
  result.insider_percentage = 100.0 * static_cast<double>(insiders)           / points.size();
  result.arc_coverage_deg   = static_cast<float>(min_arc);

  // Parameters packed as [cx1, cy1, r1,  cx2, cy2, r2,  ...]
  for (const auto& c : fitted_circles)
  {
    result.parameters.push_back(c.cx);
    result.parameters.push_back(c.cy);
    result.parameters.push_back(c.radius);
  }

  // Nodes: full outline for each circle, concatenated
  for (const auto& c : fitted_circles)
  {
    for (int deg = 0; deg <= 360; deg += 2)
    {
      double t = deg * M_PI / 180.0;
      result.nodes.push_back({ c.cx + c.radius * std::cos(t),
                             c.cy + c.radius * std::sin(t),
                             m_zmean });
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

FittingMultiCircle::CircleParams FittingMultiCircle::fit_circle_ransac(
    const std::vector<Vec3>& points,
    const std::vector<int>&  indices,
    double                   tolerance) const
{
  CircleParams best;
  best.valid = false;

  int n = static_cast<int>(indices.size());
  if (n < 3) return best;

  std::uniform_int_distribution<int> dist(0, n - 1);

  int max_inliers         = 0;
  int early_exit_threshold = static_cast<int>(m_early_exit_ratio * n);

  for (int iter = 0; iter < m_max_iterations; ++iter)
  {
    // Pick 3 unique local indices (into `indices`, not into `points`)
    int li1 = dist(m_rng);
    int li2, li3;
    do { li2 = dist(m_rng); } while (li2 == li1);
    do { li3 = dist(m_rng); } while (li3 == li1 || li3 == li2);

    const Vec3& p1 = points[indices[li1]];
    const Vec3& p2 = points[indices[li2]];
    const Vec3& p3 = points[indices[li3]];

    CircleParams circle = fit_circle_on_3_points(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
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
      if (max_inliers >= early_exit_threshold) break;
    }
  }

  if (max_inliers > 0) best.valid = true;
  return best;
}

FittingMultiCircle::CircleParams FittingMultiCircle::fit_circle_on_3_points(
    double x1, double y1,
    double x2, double y2,
    double x3, double y3) const
{
  CircleParams res;
  res.valid = false;

  double A = 2.0 * (x2 - x1);
  double B = 2.0 * (y2 - y1);
  double C = x2*x2 + y2*y2 - x1*x1 - y1*y1;
  double D = 2.0 * (x3 - x1);
  double E = 2.0 * (y3 - y1);
  double G = x3*x3 + y3*y3 - x1*x1 - y1*y1;

  double denom = A * E - B * D;
  if (std::fabs(denom) < 1e-12) return res;

  double cx = (C * E - B * G) / denom;
  double cy = (A * G - C * D) / denom;
  double r  = std::sqrt((x1 - cx) * (x1 - cx) + (y1 - cy) * (y1 - cy));

  res.cx = cx; res.cy = cy; res.radius = r; res.valid = true;
  return res;
}

double FittingMultiCircle::point_to_circle_distance(double px, double py, const CircleParams& circle) const
{
  double dx = px - circle.cx;
  double dy = py - circle.cy;
  return std::abs(std::sqrt(dx * dx + dy * dy) - circle.radius);
}

std::vector<int> FittingMultiCircle::find_inliers(
    const std::vector<Vec3>& points,
    const std::vector<int>&  candidates,
    const CircleParams&      circle,
    double                   tolerance) const
{
  std::vector<int> inliers;
  for (int idx : candidates)
    if (point_to_circle_distance(points[idx].x, points[idx].y, circle) <= tolerance)
      inliers.push_back(idx);
    return inliers;
}

} // namespace arbor::utils::fitting

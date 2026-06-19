/**
 * @file fitting_circle.cpp
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

#include "fitting_circle.h"
#include <cmath>
#include <limits>

namespace arbor::utils::fitting {

CircleFitter::CircleFitter(int max_iterations, double early_exit_ratio, unsigned seed)
  : m_max_iterations(max_iterations),
    m_early_exit_ratio(early_exit_ratio),
    m_rng(seed)
{
}

FittingResult CircleFitter::fit(const std::vector<Vec3>& points, double tolerance)
{
  FittingResult result;
  result.shape_type = "circle";

  if (points.size() < 3) { result.success = false; return result; }

  double zsum = 0.0;
  for (const auto& v : points) zsum += v.z;
  m_zmean = zsum / static_cast<double>(points.size());

  CircleModel circle = fit_circle_ransac(points, tolerance);
  if (!circle.valid) { result.success = false; return result; }

  std::vector<int> inliers = find_inliers(points, circle, tolerance);
  if (inliers.empty()) { result.success = false; return result; }

  Vec3 center_3d = calculate_3d_center(circle);

  int interiors = 0;
  for (const auto& p : points)
  {
    const double dx = p.x - circle.cx;
    const double dy = p.y - circle.cy;
    if (std::sqrt(dx*dx + dy*dy) < circle.radius - tolerance) ++interiors;
  }

  result.success             = true;
  result.center              = center_3d;
  result.radius              = static_cast<float>(circle.radius);
  result.inlier_indices      = inliers;
  result.inlier_percentage   = 100.0f * static_cast<float>(inliers.size())  / static_cast<float>(points.size());
  result.interior_percentage = 100.0f * static_cast<float>(interiors)       / static_cast<float>(points.size());
  result.arc_coverage_deg    = static_cast<float>(calculate_arc_coverage(points, inliers, center_3d));
  result.parameters          = {circle.cx, circle.cy, circle.radius};

  const double r = circle.radius;
  for (int i = 0; i <= 360; i += 2)
  {
    const double t = i * M_PI / 180.0;
    result.contour.push_back({center_3d.x + r * std::cos(t),
                              center_3d.y + r * std::sin(t),
                              m_zmean});
  }

  return result;
}

CircleModel CircleFitter::fit_circle_ransac(const std::vector<Vec3>& points, double tolerance) const
{
  CircleModel best;
  const int n = static_cast<int>(points.size());
  if (n < 3) return best;

  std::uniform_int_distribution<int> dist(0, n - 1);
  int max_inliers          = 0;
  const int early_exit_thr = static_cast<int>(m_early_exit_ratio * n);

  for (int iter = 0; iter < m_max_iterations; ++iter)
  {
    int idx1 = dist(m_rng);
    int idx2; do { idx2 = dist(m_rng); } while (idx2 == idx1);
    int idx3; do { idx3 = dist(m_rng); } while (idx3 == idx1 || idx3 == idx2);

    CircleModel circle = fit_circle_on_3_points(
      points[idx1].x, points[idx1].y,
      points[idx2].x, points[idx2].y,
      points[idx3].x, points[idx3].y);

    if (!circle.valid || circle.radius <= 0.0 || std::isnan(circle.radius)) continue;

    int inliers = 0;
    for (const auto& p : points)
      if (point_to_circle_distance(p.x, p.y, circle) < tolerance) ++inliers;

    if (inliers > max_inliers)
    {
      max_inliers = inliers;
      best = circle;
      if (max_inliers >= early_exit_thr) break;
    }
  }

  if (max_inliers > 0) best.valid = true;
  return best;
}

std::vector<int> CircleFitter::find_inliers(
    const std::vector<Vec3>& points, const CircleModel& circle, double tolerance) const
{
  std::vector<int> inliers;
  for (size_t i = 0; i < points.size(); ++i)
    if (point_to_circle_distance(points[i].x, points[i].y, circle) <= tolerance)
      inliers.push_back(static_cast<int>(i));
  return inliers;
}

Vec3 CircleFitter::calculate_3d_center(const CircleModel& circle) const
{
  return {circle.cx, circle.cy, m_zmean};
}

} // namespace arbor::utils::fitting

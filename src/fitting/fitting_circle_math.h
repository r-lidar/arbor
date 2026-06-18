/**
 * @file fitting_circle_math.h
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

#include "fitting_types.h"
#include <algorithm>
#include <cmath>

namespace arbor::utils::fitting {

/**
 * @brief Minimal representation of a fitted 2-D circle.
 */
struct CircleModel
{
  double cx     = 0.0;
  double cy     = 0.0;
  double radius = 0.0;
  bool   valid  = false;
};

/**
 * @brief Compute the unique circle passing through three 2-D points.
 *
 * Returns an invalid CircleModel when the three points are collinear
 * (determinant below 1e-12).
 */
inline CircleModel fit_circle_on_3_points(
    double x1, double y1,
    double x2, double y2,
    double x3, double y3)
{
  CircleModel res;

  const double A = 2.0 * (x2 - x1);
  const double B = 2.0 * (y2 - y1);
  const double C = x2*x2 + y2*y2 - x1*x1 - y1*y1;
  const double D = 2.0 * (x3 - x1);
  const double E = 2.0 * (y3 - y1);
  const double G = x3*x3 + y3*y3 - x1*x1 - y1*y1;

  const double denom = A * E - B * D;
  if (std::fabs(denom) < 1e-12) return res;   // collinear — no circle

  res.cx     = (C * E - B * G) / denom;
  res.cy     = (A * G - C * D) / denom;
  res.radius = std::sqrt((x1 - res.cx)*(x1 - res.cx) + (y1 - res.cy)*(y1 - res.cy));
  res.valid  = true;
  return res;
}


inline double point_to_circle_distance(double px, double py, const CircleModel& c)
{
  const double dx = px - c.cx;
  const double dy = py - c.cy;
  return std::abs(std::sqrt(dx*dx + dy*dy) - c.radius);
}

/**
 * @brief Area of the lens-shaped intersection of two circles.
 *
 * @param r1 radius of the first circle
 * @param r2 radius of the second circle
 * @param d  distance between the two centers
 *
 * Returns 0 when the circles do not overlap, and the area of the smaller
 * circle when it is fully contained inside the larger one.
 */
inline double circle_intersection_area(double r1, double r2, double d)
{
  if (r1 <= 0.0 || r2 <= 0.0) return 0.0;
  if (d >= r1 + r2) return 0.0;                            // disjoint
  if (d <= std::fabs(r1 - r2))                             // one fully inside the other
    return M_PI * std::min(r1, r2) * std::min(r1, r2);

  const double r1sq = r1 * r1;
  const double r2sq = r2 * r2;

  // Distance from circle 1's center to the radical line, and the remainder
  // to circle 2's center, then sum the two circular-segment areas.
  const double d1 = (d*d - r2sq + r1sq) / (2.0 * d);
  const double d2 = d - d1;

  // Clamp: numerical noise can push d1/r1 or d2/r2 just outside [-1, 1].
  const double a1 = std::acos(std::clamp(d1 / r1, -1.0, 1.0));
  const double a2 = std::acos(std::clamp(d2 / r2, -1.0, 1.0));

  const double area1 = r1sq * a1 - d1 * std::sqrt(std::max(0.0, r1sq - d1*d1));
  const double area2 = r2sq * a2 - d2 * std::sqrt(std::max(0.0, r2sq - d2*d2));

  return area1 + area2;
}

/**
 * @brief Overlap between two circles, relative to the area of the smaller one.
 *
 * 0.0 = the circles are disjoint.
 * 1.0 = the smaller circle is fully engulfed by the larger one.
 */
inline double circle_overlap_ratio(const CircleModel& a, const CircleModel& b)
{
  const double smaller_area = M_PI * std::min(a.radius, b.radius) * std::min(a.radius, b.radius);
  if (smaller_area <= 0.0) return 0.0;

  const double dx = a.cx - b.cx;
  const double dy = a.cy - b.cy;
  const double d  = std::sqrt(dx*dx + dy*dy);

  return circle_intersection_area(a.radius, b.radius, d) / smaller_area;
}

} // namespace arbor::utils::fitting

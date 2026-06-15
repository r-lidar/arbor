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

/**
 * @brief Unsigned distance from a 2-D point to the perimeter of a circle.
 */
inline double point_to_circle_distance(double px, double py, const CircleModel& c)
{
  const double dx = px - c.cx;
  const double dy = py - c.cy;
  return std::abs(std::sqrt(dx*dx + dy*dy) - c.radius);
}

} // namespace arbor::utils::fitting

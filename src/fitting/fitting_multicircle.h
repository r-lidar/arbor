/**
 * @file fitting_multicircle.h
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

#include "fitting_strategy.h"
#include "fitting_circle_math.h"
#include <random>
#include <vector>

namespace arbor::utils::fitting {

/**
 * @brief Fits exactly n_circles circles sequentially via RANSAC.
 *
 * Each circle is fitted on the points not already claimed by a previous
 * circle.  The combined FittingResult has:
 *   center      — area-weighted centroid of all circle centres
 *   radius      — equivalent cross-section radius: √(Σ rᵢ²)
 *   parameters  — packed as [cx₁, cy₁, r₁,  cx₂, cy₂, r₂, …]
 *   contour     — full outline of every circle, concatenated
 *   arc_coverage_deg — minimum across all circles (most conservative)
 *
 * Suitable for multi-stem cross-sections (double stem, triple stem).
 */
class MultiCircleFitter : public ShapeFitter
{
public:
  MultiCircleFitter(int n_circles = 2, int max_iterations = 1000, double early_exit_ratio = 0.9, unsigned seed = 64);
  FittingResult fit(const std::vector<Vec3>& points, double tolerance) override;

private:
  CircleModel      fit_circle_ransac(const std::vector<Vec3>& points, const std::vector<int>& indices, double tolerance) const;
  std::vector<int> find_inliers(const std::vector<Vec3>& points, const std::vector<int>& candidates, const CircleModel& circle, double tolerance) const;

  int                  m_n_circles;
  int                  m_max_iterations;
  double               m_early_exit_ratio;
  mutable std::mt19937 m_rng;
};

} // namespace arbor::utils::fitting

/**
 * @file fitting_circle.h
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
 * @brief Fits a single circle to a 2-D point cloud via RANSAC.
 */
class CircleFitter : public ShapeFitter
{
public:
  CircleFitter(int max_iterations = 1000, double early_exit_ratio = 0.9, unsigned seed = 64);
  FittingResult fit(const std::vector<Vec3>& points, double tolerance) override;

private:
  CircleModel      fit_circle_ransac(const std::vector<Vec3>& points, double tolerance) const;
  std::vector<int> find_inliers(const std::vector<Vec3>& points, const CircleModel& circle, double tolerance) const;
  Vec3             calculate_3d_center(const CircleModel& circle) const;

  int                  m_max_iterations;
  double               m_early_exit_ratio;
  mutable std::mt19937 m_rng;
};

} // namespace arbor::utils::fitting

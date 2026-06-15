/**
 * @file fitting_ellipse.h
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
#include <random>
#include <vector>

namespace arbor::utils::fitting {

/**
 * @brief Fits an ellipse to a 2-D point cloud via algebraic least-squares with RANSAC.
 *
 * The algebraic form ax² + bxy + cy² + dx + ey + f = 0 is solved by finding the
 * eigenvector of the 6×6 scatter matrix that corresponds to the smallest eigenvalue
 * (Jacobi decomposition).  RANSAC is used for robustness against outliers.
 */
class EllipseFitter : public ShapeFitter
{
public:
  EllipseFitter(int max_iterations = 1000, int min_inliers = 10, unsigned seed = 42);
  FittingResult fit(const std::vector<Vec3>& points, double tolerance) override;

private:
  struct EllipseParams
  {
    double a, b, c, d, e, f;
    bool valid = false;
  };

  struct EllipseGeometry
  {
    double cx, cy, major, minor, angle;
    bool valid = false;
  };

  EllipseParams       fit_ellipse_algebraic(const std::vector<Vec3>& pts) const;
  std::vector<double> calculate_distances(const std::vector<Vec3>& pts, const EllipseParams& params) const;
  EllipseGeometry     get_ellipse_geometry(const EllipseParams& params) const;

  int                  m_max_iterations;
  int                  m_min_inliers;
  mutable std::mt19937 m_rng;
};

} // namespace arbor::utils::fitting

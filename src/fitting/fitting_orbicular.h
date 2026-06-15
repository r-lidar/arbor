/**
 * @file fitting_orbicular.h
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
#include <array>
#include <memory>
#include <vector>

namespace arbor::utils::fitting {

/**
 * @brief Public-facing orchestrator for tree-stem cross-section fitting.
 *
 * Accepts a 3-D point cloud that lies roughly on a planar cross-section of a
 * tree stem or branch.  Internally the cloud is translated to the local origin
 * and rotated so the cylinder axis aligns with Z, reducing the problem to 2-D.
 * One or more ShapeFitter strategies are then tried; the one with the highest
 * inlier percentage is returned, with its center and contour transformed back
 * to the original coordinate frame.
 *
 * Typical usage:
 * @code
 *   CrossSectionFitter fitter;
 *   fitter.set_axis(bottom, top);
 *   for (auto& p : cloud) fitter.add_point(p.x, p.y, p.z);
 *   auto result = fitter.fit(0.03, FitMode::Standard);
 * @endcode
 */
class CrossSectionFitter
{
public:
  CrossSectionFitter();

  /** Define the cylinder axis.  Points added afterward are projected relative to it. */
  void set_axis(const Vec3& from, const Vec3& to);

  void add_point(double x, double y, double z);
  void clear();

  FittingResult fit(double tolerance = 0.03, FitMode flags = FitMode::Basic);

private:
  Vec3 m_origin;
  std::vector<Vec3> m_points;
  std::array<std::array<double, 3>, 3> R;
  std::array<std::array<double, 3>, 3> R_inv;

  void compute_rotation_matrix(const Vec3& from, const Vec3& to);
  void apply_rotation(Vec3& p);
  void reverse_rotation(Vec3& p);
};

} // namespace arbor::utils::fitting

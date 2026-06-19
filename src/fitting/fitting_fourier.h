/**
 * @file fitting_fourier.h
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
#include <vector>

namespace arbor::utils::fitting {

/**
 * @brief Fits an irregular closed contour to a 2-D point cloud using a
 *        truncated Fourier series in polar coordinates.
 *
 * The contour is expressed as:
 *   r(θ) = a₀ + Σₖ [ aₖ cos(kθ) + bₖ sin(kθ) ],  k = 1 … K
 *
 * The number of harmonics K is set at construction time.
 * Gaps in angular coverage are filled with the median radius before fitting
 * so that missing sectors do not bias the coefficients.
 *
 * Suitable for complex cross-sections: leaning stems, irregular buttresses.
 */
class FourierFitter : public ShapeFitter
{
public:
  FourierFitter(int K = 5, double step_deg = 10.0);
  FittingResult fit(const std::vector<Vec3>& points, double tolerance) override;

private:
  int    m_K;
  double m_step_deg;

  Vec3        calculate_centroid(const std::vector<Vec3>& points) const;
  PolarData   to_polar(const std::vector<Vec3>& points, const Vec3& center) const;
  double      calculate_median(std::vector<double> values) const;
  PolarData   inject_missing_angles(const PolarData& polar, double fill_radius) const;

  std::vector<double> build_fourier_matrix(const std::vector<double>& theta) const;
  std::vector<double> fit_fourier(const std::vector<double>& theta, const std::vector<double>& r) const;
  double              evaluate_fourier(double theta, const std::vector<double>& coefficients) const;
};

} // namespace arbor::utils::fitting

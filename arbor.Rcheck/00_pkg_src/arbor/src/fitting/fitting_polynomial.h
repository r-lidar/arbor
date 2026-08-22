/**
 * @file fitting_polynomial.h
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

#include <vector>
#include <utility>
#include <limits>
#include <cmath>
#include <stdexcept>

namespace arbor::fitting {

// Fits a quadratic  R(x) = intersect + a*x + b*x²
// to (subtree_length, radius - intersect) data points.
//
// In DecreasingOnly mode the derivative dR/dx = a + 2*b*x is constrained
// to be <= 0 over the entire data domain, ensuring radii never increase
// from tip toward root.
// PolynomialFitting.h
class PolynomialFitting
{
public:
  enum class Mode { Unconstrained, DecreasingOnly };
  PolynomialFitting(const std::vector<std::pair<float, float>>& pts, float intersect, Mode mode = Mode::DecreasingOnly);
  float predict(float x) const;
  bool valid = false;

private:
  float intersect;
  float a = 0.0;
  float b = 0.0;

  float x_min =  std::numeric_limits<float>::max();
  float x_max = -std::numeric_limits<float>::max();
  double sum_x2 = 0, sum_x3 = 0, sum_x4 = 0, sum_xy = 0, sum_x2y = 0;

  void accumulate(const std::vector<std::pair<float, float>>& pts);
  bool solve_unconstrained();
  bool solve_active(float x_c);
  bool is_feasible() const;
  float rss(const std::vector<std::pair<float, float>>& pts, float a_, float b_) const;
};

}

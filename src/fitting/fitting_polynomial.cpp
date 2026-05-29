/**
 * @file fitting_polynomial.cpp
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

#include "fitting_polynomial.h"

namespace arbor::fitting {

PolynomialFitting::PolynomialFitting(const std::vector<std::pair<float, float>>& pts, float intersect_, Mode mode) : intersect(intersect_)
{
  if (pts.size() < 3) return;

  accumulate(pts);

  if (mode == Mode::Unconstrained)
  {
    valid = solve_unconstrained();
    return;
  }

  // DecreasingOnly: enumerate KKT cases, keep best feasible RSS
  // Case 3 fallback: a=0, b=0 => dR/dx=0, always satisfies <=0
  float best = rss(pts, 0.0f, 0.0f);
  valid = true;

  auto try_candidate = [&](float a_, float b_)
  {
    float r = rss(pts, a_, b_);
    if (r < best) { best = r; a = a_; b = b_; }
  };

  float ca, cb;

  // Case 0: unconstrained
  ca = a; cb = b;
  if (solve_unconstrained() && is_feasible())
    try_candidate(a, b);

  // Case 1: active at x_min
  ca = a; cb = b;
  if (solve_active(x_min) && is_feasible())
    try_candidate(a, b);
  a = ca; b = cb;

  // Case 2: active at x_max
  if (solve_active(x_max) && is_feasible())
    try_candidate(a, b);
}

float PolynomialFitting::predict(float x) const
{
  float pred = intersect + a * x + b * x * x;
  return pred < 0.0 ? 0.0 : pred;
}

void PolynomialFitting::accumulate(const std::vector<std::pair<float, float>>& pts)
{
  for (const auto& [x, y] : pts)
  {
    float x2 = x*x, x3 = x2*x, x4 = x3*x;
    sum_x2  += x2;
    sum_x3  += x3;
    sum_x4  += x4;
    sum_xy  += x  * y;
    sum_x2y += x2 * y;
    if (x < x_min) x_min = x;
    if (x > x_max) x_max = x;
  }
}

bool PolynomialFitting::solve_unconstrained()
{
  double det = sum_x2 * sum_x4 - sum_x3 * sum_x3;
  if (std::abs(det) < 1e-12) return false;
  a = (sum_xy * sum_x4  - sum_x2y * sum_x3) / det;
  b = (sum_x2 * sum_x2y - sum_xy  * sum_x3) / det;
  return true;
}

bool PolynomialFitting::solve_active(float x_c)
{
  // Substitute a = -2*b*x_c, let z = x² - 2*x_c*x
  // => b = Σ(z*y) / Σ(z²), derived from stored sums
  double sz2 = sum_x4 - 4.0*x_c*sum_x3 + 4.0*x_c*x_c*sum_x2;
  double szy = sum_x2y - 2.0*x_c*sum_xy;
  if (std::abs(sz2) < 1e-12) return false;
  b = szy / sz2;
  a = -2.0 * b * x_c;
  return true;
}

bool PolynomialFitting::is_feasible() const
{
  // dR/dx = a + 2*b*x must be <= 0 at both endpoints
  const double eps = 1e-9;
  return (a + 2.0*b*x_min <= eps) && (a + 2.0*b*x_max <= eps);
}

float PolynomialFitting::rss(const std::vector<std::pair<float, float>>& pts, float a_, float b_) const
{
  double s = 0.0;
  for (const auto& [x, y] : pts)
  {
    double r = a_*x + b_*x*x - y;
    s += r*r;
  }
  return s;
}

}

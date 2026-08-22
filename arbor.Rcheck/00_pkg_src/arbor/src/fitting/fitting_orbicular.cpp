/**
 * @file fitting_orbicular.cpp
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

#include "fitting_orbicular.h"
#include "fitting_circle.h"
#include "fitting_ellipse.h"
#include "fitting_fourier.h"
#include "fitting_multicircle.h"
#include <algorithm>
#include <cmath>
#include <set>

namespace arbor::utils::fitting {

double ShapeFitter::calculate_arc_coverage(
    const std::vector<Vec3>& points,
    const std::vector<int>&  inliers,
    const Vec3&              center)
{
  if (inliers.empty()) return 0.0;

  const double bin_size = 10.0;
  const int total_bins = 36; // 360 / 10

  // Use a vector of bools as a bitmask for occupied bins (0 to 35)
  std::vector<bool> has_point(total_bins, false);
  int unique_bin_count = 0;

  for (int idx : inliers)
  {
    const auto& p = points[idx];
    double angle = std::atan2(p.y - center.y, p.x - center.x) * 180.0 / M_PI;
    if (angle < 0.0) angle += 360.0;

    // Integer division to map 0.0-359.999 to 0-35 indices safely
    int bin_idx = static_cast<int>(std::lround(angle / bin_size)) % total_bins;
    if (bin_idx >= total_bins) bin_idx = total_bins - 1; // Guard against rounding edge cases

    if (!has_point[bin_idx]) {
      has_point[bin_idx] = true;
      ++unique_bin_count;
    }
  }

  if (unique_bin_count == 0) return 0.0;
  if (unique_bin_count == total_bins) return 360.0;

  // Find the longest consecutive chain of occupied bins
  int max_chain = 0;
  int current_chain = 0;

  // Run the loop up to 2 * total_bins to seamlessly handle the 360->0 wraparound
  for (int i = 0; i < 2 * total_bins; ++i)
  {
    if (has_point[i % total_bins])
    {
      ++current_chain;
      if (current_chain > max_chain) {
        max_chain = current_chain;
      }
    }
    else
    {
      current_chain = 0;
    }
  }

  // Cap the max chain at total_bins just in case
  if (max_chain > total_bins) max_chain = total_bins;

  return max_chain * bin_size;
}

CrossSectionFitter::CrossSectionFitter()
{
  R = {{
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}
  }};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      R_inv[i][j] = R[j][i];
}

void CrossSectionFitter::set_axis(const Vec3& from, const Vec3& to)
{
  m_origin = from;
  const Vec3 axis_dir = (to - from).normalized();
  const Vec3 z_axis(0.0, 0.0, 1.0);
  compute_rotation_matrix(axis_dir, z_axis);
}

void CrossSectionFitter::add_point(double x, double y, double z)
{
  Vec3 p = {x - m_origin.x, y - m_origin.y, z - m_origin.z};
  apply_rotation(p);
  m_points.push_back(p);
}

void CrossSectionFitter::clear() { m_points.clear(); }

FittingResult CrossSectionFitter::fit(double tolerance, FitMode flags)
{
  FittingResult best_result;
  std::vector<std::unique_ptr<ShapeFitter>> strategies;

  if (has(flags, FitMode::Circle))       strategies.push_back(std::make_unique<CircleFitter>());
  if (has(flags, FitMode::Ellipse))      strategies.push_back(std::make_unique<EllipseFitter>());
  if (has(flags, FitMode::Fourier5))     strategies.push_back(std::make_unique<FourierFitter>(5));
  if (has(flags, FitMode::MultiCircle2)) strategies.push_back(std::make_unique<MultiCircleFitter>(2));
  if (has(flags, FitMode::MultiCircle3)) strategies.push_back(std::make_unique<MultiCircleFitter>(3));
  if (has(flags, FitMode::Fourier10))    strategies.push_back(std::make_unique<FourierFitter>(10));

  for (auto& strategy : strategies)
  {
    FittingResult current = strategy->fit(m_points, tolerance);
    if (current.success && current.inlier_percentage > best_result.inlier_percentage * 1.1)
      best_result = std::move(current);
  }

  // Transform the best result back to the original coordinate frame.
  reverse_rotation(best_result.center);
  best_result.center.x += m_origin.x;
  best_result.center.y += m_origin.y;
  best_result.center.z += m_origin.z;

  for (auto& v : best_result.contour)
  {
    reverse_rotation(v);
    v.x += m_origin.x;
    v.y += m_origin.y;
    v.z += m_origin.z;
  }

  return best_result;
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

void CrossSectionFitter::compute_rotation_matrix(const Vec3& from, const Vec3& to)
{
  const Vec3 a = from.normalized();
  const Vec3 b = to.normalized();

  Vec3 v  = a.cross(b);     // rotation axis scaled by sin(angle)
  double s = v.length();    // sin(angle)
  double c = a.dot(b);      // cos(angle)

  if (s < 1e-9)
  {
    if (c > 0) return;       // already aligned — R stays identity

    // Vectors are 180° apart: pick an arbitrary perpendicular axis
    const Vec3 axis = (std::fabs(a.x) > 0.9) ? Vec3(0, 1, 0) : Vec3(1, 0, 0);
    v = a.cross(axis).normalized();
    s = 0.0;
    c = -1.0;
  }
  else
  {
    v = v / s;               // normalize rotation axis
  }

  // Rodrigues' rotation formula: R = I + sin(θ)·K + (1 − cos(θ))·K²
  const double vx = v.x, vy = v.y, vz = v.z;
  const double K[3][3] = {
    { 0,  -vz,  vy},
    { vz,  0,  -vx},
    {-vy,  vx,  0 }
  };

  const double one_minus_c = 1.0 - c;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      const double identity  = (i == j) ? 1.0 : 0.0;
      const double k_squared = K[i][0]*K[0][j] + K[i][1]*K[1][j] + K[i][2]*K[2][j];
      R[i][j] = identity + K[i][j]*s + k_squared*one_minus_c;
    }
  }

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      R_inv[i][j] = R[j][i];   // R is orthogonal: R⁻¹ = Rᵀ
}

void CrossSectionFitter::apply_rotation(Vec3& p)
{
  const double x = R[0][0]*p.x + R[0][1]*p.y + R[0][2]*p.z;
  const double y = R[1][0]*p.x + R[1][1]*p.y + R[1][2]*p.z;
  const double z = R[2][0]*p.x + R[2][1]*p.y + R[2][2]*p.z;
  p.x = x; p.y = y; p.z = z;
}

void CrossSectionFitter::reverse_rotation(Vec3& p)
{
  const double x = R_inv[0][0]*p.x + R_inv[0][1]*p.y + R_inv[0][2]*p.z;
  const double y = R_inv[1][0]*p.x + R_inv[1][1]*p.y + R_inv[1][2]*p.z;
  const double z = R_inv[2][0]*p.x + R_inv[2][1]*p.y + R_inv[2][2]*p.z;
  p.x = x; p.y = y; p.z = z;
}

} // namespace arbor::utils::fitting

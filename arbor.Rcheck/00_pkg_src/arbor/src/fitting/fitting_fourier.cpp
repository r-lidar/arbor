/**
 * @file fitting_fourier.cpp
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

#include "fitting_fourier.h"
#include <algorithm>
#include <cmath>
#include <numeric>

// ---------------------------------------------------------------------------
// Minimal dense linear-algebra helpers (row-major: M[i][j] = M[i*cols + j])
// ---------------------------------------------------------------------------
namespace {

using Mat = std::vector<double>;
using Vec = std::vector<double>;

Mat mat_zeros(int rows, int cols) { return Mat(rows * cols, 0.0); }

// C = Aᵀ A  (A is m×n)
Mat mat_ata(const Mat& A, int m, int n)
{
  Mat C = mat_zeros(n, n);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k)
        C[j*n + k] += A[i*n + j] * A[i*n + k];
  return C;
}

// b = Aᵀ y  (A is m×n, y is length m)
Vec mat_aty(const Mat& A, const Vec& y, int m, int n)
{
  Vec b(n, 0.0);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      b[j] += A[i*n + j] * y[i];
  return b;
}

// Solve (S + eps·I) x = b for symmetric n×n S via Cholesky (LL^T).
// eps provides Tikhonov regularisation against near-singular systems.
Vec solve_symmetric(Mat S, Vec b, int n, double eps = 1e-9)
{
  for (int i = 0; i < n; ++i) S[i*n + i] += eps;

  Mat L = mat_zeros(n, n);
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      double s = S[i*n + j];
      for (int k = 0; k < j; ++k) s -= L[i*n + k] * L[j*n + k];
      if (i == j) { L[i*n + i] = std::sqrt(std::max(s, 0.0)); }
      else        { L[i*n + j] = (L[j*n + j] > 0.0) ? s / L[j*n + j] : 0.0; }
    }
  }

  // Forward substitution: L z = b
  Vec z(n, 0.0);
  for (int i = 0; i < n; ++i)
  {
    double s = b[i];
    for (int k = 0; k < i; ++k) s -= L[i*n + k] * z[k];
    z[i] = (L[i*n + i] > 0.0) ? s / L[i*n + i] : 0.0;
  }

  // Back substitution: Lᵀ x = z
  Vec x(n, 0.0);
  for (int i = n - 1; i >= 0; --i)
  {
    double s = z[i];
    for (int k = i + 1; k < n; ++k) s -= L[k*n + i] * x[k];
    x[i] = (L[i*n + i] > 0.0) ? s / L[i*n + i] : 0.0;
  }

  return x;
}

} // anonymous namespace

namespace arbor::utils::fitting {

FourierFitter::FourierFitter(int K, double step_deg) : m_K(K), m_step_deg(step_deg)
{
}

FittingResult FourierFitter::fit(const std::vector<Vec3>& points, double tolerance)
{
  FittingResult result;
  result.shape_type = "fourier" + std::to_string(m_K);
  if (points.size() < 3) return result;

  double zsum = 0.0;
  for (const auto& v : points) zsum += v.z;
  m_zmean = zsum / static_cast<double>(points.size());

  result.center = calculate_centroid(points);

  const PolarData polar     = to_polar(points, result.center);
  const double fill_radius  = calculate_median(polar.r);
  const PolarData augmented = inject_missing_angles(polar, fill_radius);

  const std::vector<double> coefficients = fit_fourier(augmented.theta, augmented.r);

  for (size_t i = 0; i < polar.theta.size(); ++i)
  {
    if (std::abs(polar.r[i] - evaluate_fourier(polar.theta[i], coefficients)) < tolerance)
    {
      result.inlier_indices.push_back(static_cast<int>(i));
    }
  }

  result.success            = !result.inlier_indices.empty();
  result.inlier_percentage  = 100.0f * static_cast<float>(result.inlier_indices.size()) / static_cast<float>(points.size());
  result.arc_coverage_deg   = static_cast<float>(calculate_arc_coverage(points, result.inlier_indices, result.center));
  if (result.arc_coverage_deg < 300.0f) result.success = false;

  for (double coeff : coefficients) result.parameters.push_back(coeff);

  for (int i = 0; i <= 360; i += 2)
  {
    const double angle_rad = i * M_PI / 180.0;
    const double radius    = evaluate_fourier(angle_rad, coefficients);
    const double xcontour  = result.center.x + radius * std::cos(angle_rad);
    const double ycontour  = result.center.y + radius * std::sin(angle_rad);
    result.contour.push_back({xcontour, ycontour, m_zmean});
  }

  // Recompute center and radius from the actual contour polygon via the
  // shoelace / centroid-of-polygon formula.
  const size_t n     = result.contour.size();
  double A           = 0.0;
  double Cx          = 0.0;
  double Cy          = 0.0;
  const double ref_x = result.center.x;
  const double ref_y = result.center.y;

  for (size_t i = 0; i < n; ++i)
  {
    const auto& p0  = result.contour[i];
    const auto& p1  = result.contour[(i + 1) % n];
    const double x0 = p0.x - ref_x,  y0 = p0.y - ref_y;
    const double x1 = p1.x - ref_x,  y1 = p1.y - ref_y;
    const double cross = x0*y1 - x1*y0;
    A  += cross;
    Cx += (x0 + x1) * cross;
    Cy += (y0 + y1) * cross;
  }
  A *= 0.5;

  if (std::abs(A) < 1e-9)
  {
    result.radius = 0.0f;
  }
  else
  {
    result.radius = static_cast<float>(std::sqrt(std::abs(A) / M_PI));
    Cx /= (6.0 * A);
    Cy /= (6.0 * A);
    result.center = {Cx + ref_x, Cy + ref_y, m_zmean};
  }

  return result;
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

Vec3 FourierFitter::calculate_centroid(const std::vector<Vec3>& points) const
{
  double sx = 0, sy = 0, sz = 0;
  for (const auto& p : points) { sx += p.x; sy += p.y; sz += p.z; }
  const double n = static_cast<double>(points.size());
  return {sx/n, sy/n, sz/n};
}

PolarData FourierFitter::to_polar(const std::vector<Vec3>& points, const Vec3& center) const
{
  PolarData p;
  p.theta.reserve(points.size());
  p.r.reserve(points.size());
  for (const auto& pt : points)
  {
    const double dx = pt.x - center.x, dy = pt.y - center.y;
    p.theta.push_back(std::fmod(std::atan2(dy, dx) + 2.0*M_PI, 2.0*M_PI));
    p.r.push_back(std::sqrt(dx*dx + dy*dy));
  }
  return p;
}

double FourierFitter::calculate_median(std::vector<double> v) const
{
  if (v.empty()) return 0.0;
  std::nth_element(v.begin(), v.begin() + v.size()/2, v.end());
  return v[v.size()/2];
}

PolarData FourierFitter::inject_missing_angles(const PolarData& polar, double fill_radius) const
{
  const double step_rad = m_step_deg * M_PI / 180.0;

  std::vector<size_t> idx(polar.theta.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b){ return polar.theta[a] < polar.theta[b]; });

  PolarData res;
  for (size_t i = 0; i < idx.size(); ++i)
  {
    const size_t curr = idx[i];
    const size_t next = idx[(i + 1) % idx.size()];
    res.theta.push_back(polar.theta[curr]);
    res.r.push_back(polar.r[curr]);

    const bool   last = (i == idx.size() - 1);
    const double diff = last
      ? (2.0*M_PI - polar.theta[curr] + polar.theta[next])
      : (polar.theta[next] - polar.theta[curr]);

    if (diff > step_rad)
    {
      const double end = last ? 2.0*M_PI + polar.theta[next] : polar.theta[next];
      for (double a = polar.theta[curr] + step_rad; a < end; a += step_rad)
      {
        res.theta.push_back(std::fmod(a, 2.0*M_PI));
        res.r.push_back(fill_radius);
      }
    }
  }
  return res;
}

// Build the Fourier design matrix (row-major, m × (1 + 2K)):
//   col 0      : 1
//   col 2k-1   : cos(k·θ)
//   col 2k     : sin(k·θ)
std::vector<double> FourierFitter::build_fourier_matrix(const std::vector<double>& theta) const
{
  const int m = static_cast<int>(theta.size());
  const int n = 1 + 2 * m_K;
  Mat X = mat_zeros(m, n);

  for (int i = 0; i < m; ++i)
  {
    X[i*n + 0] = 1.0;
    for (int k = 1; k <= m_K; ++k)
    {
      X[i*n + (2*k - 1)] = std::cos(k * theta[i]);
      X[i*n + (2*k    )] = std::sin(k * theta[i]);
    }
  }
  return X;
}

std::vector<double> FourierFitter::fit_fourier(const std::vector<double>& theta, const std::vector<double>& r) const
{
  const int m = static_cast<int>(theta.size());
  const int n = 1 + 2 * m_K;
  const Mat X   = build_fourier_matrix(theta);
  const Mat XtX = mat_ata(X, m, n);
  const Vec Xty = mat_aty(X, r, m, n);
  return solve_symmetric(XtX, Xty, n);
}

double FourierFitter::evaluate_fourier(double theta, const std::vector<double>& coeffs) const
{
  double r = coeffs[0];
  for (int k = 1; k <= m_K; ++k)
  {
    r += coeffs[2*k - 1] * std::cos(k * theta) + coeffs[2*k] * std::sin(k * theta);
  }
  return r;
}

} // namespace arbor::utils::fitting

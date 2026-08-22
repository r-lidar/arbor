/**
 * @file fitting_ellipse.cpp
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

#include "fitting_ellipse.h"
#include <array>
#include <cmath>
#include <limits>

namespace arbor::utils::fitting {

EllipseFitter::EllipseFitter(int max_iterations, int min_inliers, unsigned seed)
  : m_max_iterations(max_iterations), m_min_inliers(min_inliers), m_rng(seed)
{
}

EllipseFitter::EllipseParams EllipseFitter::fit_ellipse_algebraic(const std::vector<Vec3>& pts) const
{
  EllipseParams res;
  if (pts.size() < 5) return res;

  double mean_x = 0.0;
  double mean_y = 0.0;
  for (const auto& p : pts)
  {
    mean_x += p.x;
    mean_y += p.y;
  }
  mean_x /= static_cast<double>(pts.size());
  mean_y /= static_cast<double>(pts.size());

  double sq_sum = 0.0;
  for (const auto& p : pts)
  {
    const double dx = p.x - mean_x;
    const double dy = p.y - mean_y;
    sq_sum += dx*dx + dy*dy;
  }

  const double scale = std::sqrt(sq_sum / static_cast<double>(pts.size()));
  if (!std::isfinite(scale) || scale <= std::numeric_limits<double>::epsilon()) return res;

  // Build the 6×6 scatter matrix S = Dᵀ D where each row of D is
  // [x², xy, y², x, y, 1].  The ellipse coefficients are the eigenvector
  // of S with the smallest eigenvalue.
  const int N = 6;
  double S[N][N] = {};

  for (const auto& p : pts)
  {
    const double x = (p.x - mean_x) / scale;
    const double y = (p.y - mean_y) / scale;
    const double row[6] = {x*x, x*y, y*y, x, y, 1.0};
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        S[i][j] += row[i] * row[j];
  }

  // Jacobi eigendecomposition of the symmetric 6×6 matrix S.
  // V accumulates the eigenvectors (columns); A converges to diagonal.
  double V[N][N] = {};
  double A[N][N] = {};
  for (int i = 0; i < N; ++i) { V[i][i] = 1.0; for (int j = 0; j < N; ++j) A[i][j] = S[i][j]; }

  for (int sweep = 0; sweep < 100; ++sweep)
  {
    double off = 0.0;
    for (int i = 0; i < N; ++i)
      for (int j = i+1; j < N; ++j)
        off += A[i][j] * A[i][j];
    if (off < 1e-24) break;

    for (int p2 = 0; p2 < N-1; ++p2)
    {
      for (int q = p2+1; q < N; ++q)
      {
        if (std::fabs(A[p2][q]) < 1e-15) continue;

        const double tau = (A[q][q] - A[p2][p2]) / (2.0 * A[p2][q]);
        const double t   = (tau >= 0 ? 1.0 : -1.0) / (std::fabs(tau) + std::sqrt(1.0 + tau*tau));
        const double c   = 1.0 / std::sqrt(1.0 + t*t);
        const double s   = t * c;

        const double App = A[p2][p2], Aqq = A[q][q], Apq = A[p2][q];
        A[p2][p2] = App - t*Apq;
        A[q][q]   = Aqq + t*Apq;
        A[p2][q]  = A[q][p2] = 0.0;

        for (int r = 0; r < N; ++r)
        {
          if (r == p2 || r == q) continue;
          const double Arp = A[r][p2], Arq = A[r][q];
          A[r][p2] = A[p2][r] = c*Arp - s*Arq;
          A[r][q]  = A[q][r]  = s*Arp + c*Arq;
        }
        for (int r = 0; r < N; ++r)
        {
          const double Vrp = V[r][p2], Vrq = V[r][q];
          V[r][p2] = c*Vrp - s*Vrq;
          V[r][q]  = s*Vrp + c*Vrq;
        }
      }
    }
  }

  // Find the column of V whose eigenvalue (diagonal of A) is smallest.
  // BUG FIX: min_col = i must be inside the if-block; previously it was
  // unconditionally executed, always yielding min_col = N-1 = 5.
  int    min_col = 0;
  double min_val = std::fabs(A[0][0]);
  for (int i = 1; i < N; ++i)
  {
    if (std::fabs(A[i][i]) < min_val)
    {
      min_val = std::fabs(A[i][i]);
      min_col = i;               // ← was outside the if-block before
    }
  }

  const double an = V[0][min_col];
  const double bn = V[1][min_col];
  const double cn = V[2][min_col];
  const double dn = V[3][min_col];
  const double en = V[4][min_col];
  const double fn = V[5][min_col];

  res.a = an;
  res.b = bn;
  res.c = cn;
  res.d = -2.0*an*mean_x - bn*mean_y + dn*scale;
  res.e = -bn*mean_x - 2.0*cn*mean_y + en*scale;
  res.f = an*mean_x*mean_x + bn*mean_x*mean_y + cn*mean_y*mean_y
        - dn*scale*mean_x - en*scale*mean_y + fn*scale*scale;
  res.valid = (res.b*res.b - 4.0*res.a*res.c < 0);
  return res;
}

std::vector<double> EllipseFitter::calculate_distances(const std::vector<Vec3>& pts, const EllipseParams& p) const
{
  std::vector<double> dists;
  dists.reserve(pts.size());
  for (const auto& pt : pts)
  {
    const double f  = p.a*pt.x*pt.x + p.b*pt.x*pt.y + p.c*pt.y*pt.y + p.d*pt.x + p.e*pt.y + p.f;
    const double gx = 2.0*p.a*pt.x + p.b*pt.y + p.d;
    const double gy = p.b*pt.x + 2.0*p.c*pt.y + p.e;
    dists.push_back(std::abs(f) / std::sqrt(gx*gx + gy*gy + 1e-12));
  }
  return dists;
}

EllipseFitter::EllipseGeometry EllipseFitter::get_ellipse_geometry(const EllipseParams& p) const
{
  EllipseGeometry g;
  if (!std::isfinite(p.a) || !std::isfinite(p.b) || !std::isfinite(p.c)
      || !std::isfinite(p.d) || !std::isfinite(p.e) || !std::isfinite(p.f))
    return g;

  const double det = p.b*p.b - 4.0*p.a*p.c;
  if (det >= 0) return g;

  g.cx    = (2.0*p.c*p.d - p.b*p.e) / det;
  g.cy    = (2.0*p.a*p.e - p.b*p.d) / det;
  g.angle = 0.5 * std::atan2(p.b, p.a - p.c);

  const double up   = 2.0 * (p.a*p.e*p.e + p.c*p.d*p.d + p.f*p.b*p.b - p.b*p.d*p.e - 4.0*p.a*p.c*p.f);
  const double root = std::sqrt(std::pow(p.a - p.c, 2) + p.b*p.b);
  const double major_sq = up / (det * (p.a + p.c - root));
  const double minor_sq = up / (det * (p.a + p.c + root));
  if (!std::isfinite(g.cx) || !std::isfinite(g.cy) || !std::isfinite(g.angle)
      || !std::isfinite(major_sq) || !std::isfinite(minor_sq)
      || major_sq <= 0.0 || minor_sq <= 0.0)
    return g;

  g.major = std::sqrt(major_sq);
  g.minor = std::sqrt(minor_sq);
  if (g.minor > g.major) std::swap(g.minor, g.major);
  if (!std::isfinite(g.major) || !std::isfinite(g.minor) || g.major <= 0.0 || g.minor <= 0.0)
    return g;

  g.valid = true;
  return g;
}

FittingResult EllipseFitter::fit(const std::vector<Vec3>& points, double tolerance)
{
  FittingResult best;
  best.shape_type = "ellipse";
  if (points.size() < 5) return best;

  double zsum = 0.0;
  for (const auto& v : points) zsum += v.z;
  m_zmean = zsum / static_cast<double>(points.size());

  //std::uniform_int_distribution<size_t> dist(0, points.size() - 1);
  deterministic_uniform_int_distribution dist(0, points.size() - 1);
  for (int i = 0; i < m_max_iterations; ++i)
  {
    std::vector<Vec3> sample;
    sample.reserve(5);
    std::array<size_t, 5> sample_indices{};
    size_t unique_count = 0;
    while (unique_count < sample_indices.size())
    {
      const size_t candidate = dist(m_rng);
      bool duplicate = false;
      for (size_t j = 0; j < unique_count; ++j)
      {
        if (sample_indices[j] == candidate)
        {
          duplicate = true;
          break;
        }
      }
      if (duplicate) continue;
      sample_indices[unique_count++] = candidate;
      sample.push_back(points[candidate]);
    }

    const auto p = fit_ellipse_algebraic(sample);
    if (!p.valid) continue;

    const auto d = calculate_distances(points, p);
    std::vector<int> inliers;
    for (size_t j = 0; j < d.size(); ++j)
      if (d[j] < tolerance) inliers.push_back(static_cast<int>(j));

    const auto g = get_ellipse_geometry(p);
    if (g.valid
        && inliers.size() >= static_cast<size_t>(m_min_inliers)
        && inliers.size() > best.inlier_indices.size())
    {
      best.inlier_indices = inliers;
      best.success        = g.valid && (g.major < 3.0 * g.minor);
    }
  }

  if (best.success)
  {
    std::vector<Vec3> inlier_pts;
    inlier_pts.reserve(best.inlier_indices.size());
    for (int idx : best.inlier_indices) inlier_pts.push_back(points[idx]);

    const auto final_p = fit_ellipse_algebraic(inlier_pts);
    const auto g       = get_ellipse_geometry(final_p);

    if (!final_p.valid || !g.valid)
    {
      best.success = false;
      return best;
    }

    const auto final_d = calculate_distances(points, final_p);
    std::vector<int> final_inliers;
    final_inliers.reserve(points.size());
    for (size_t i = 0; i < final_d.size(); ++i)
      if (final_d[i] < tolerance) final_inliers.push_back(static_cast<int>(i));

    if (final_inliers.size() < static_cast<size_t>(m_min_inliers))
    {
      best.success = false;
      return best;
    }

    best.radius           = static_cast<float>(std::sqrt(g.major * g.minor));
    best.center           = {g.cx, g.cy, m_zmean};
    best.parameters       = {final_p.a, final_p.b, final_p.c, final_p.d, final_p.e, final_p.f,
                              g.major, g.minor, g.angle};
    best.inlier_indices    = std::move(final_inliers);
    best.inlier_percentage = 100.0f * static_cast<float>(best.inlier_indices.size()) / static_cast<float>(points.size());
    best.arc_coverage_deg  = static_cast<float>(calculate_arc_coverage(points, best.inlier_indices, best.center));

    const double a     = g.major;
    const double b     = g.minor;
    const double theta = g.angle + M_PI / 2.0;
    const double cos_t = std::cos(theta);
    const double sin_t = std::sin(theta);

    for (int i = 0; i <= 360; i += 2)
    {
      const double t  = i * M_PI / 180.0;
      const double xr = a * std::cos(t);
      const double yr = b * std::sin(t);
      best.contour.push_back({g.cx + xr*cos_t - yr*sin_t,
                              g.cy + xr*sin_t + yr*cos_t,
                              m_zmean});
    }
  }

  return best;
}

} // namespace arbor::utils::fitting

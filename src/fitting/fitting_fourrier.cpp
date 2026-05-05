#include "fitting.h"
#include <numeric>
#include <cmath>
#include <stdexcept>
#include <algorithm>

// ---------------------------------------------------------------------------
// Minimal dense linear-algebra helpers (replace Eigen)
// All matrices are row-major: M[i][j] = M[i * cols + j]
// ---------------------------------------------------------------------------
namespace {

using Mat = std::vector<double>;  // flat row-major
  using Vec = std::vector<double>;

  Mat mat_zeros(int rows, int cols)
  {
    return Mat(rows * cols, 0.0);
  }

  // C = A^T * A   (A is m x n)
  Mat mat_ata(const Mat& A, int m, int n)
  {
    Mat C = mat_zeros(n, n);
    for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j)
        for (int k = 0; k < n; ++k)
          C[j * n + k] += A[i * n + j] * A[i * n + k];
    return C;
  }

  // b = A^T * y   (A is m x n, y is m)
  Vec mat_aty(const Mat& A, const Vec& y, int m, int n)
  {
    Vec b(n, 0.0);
    for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j)
        b[j] += A[i * n + j] * y[i];
    return b;
  }

  // Solve (S + eps*I) x = b for a symmetric n x n matrix S using
  // Cholesky decomposition (LL^T).  eps provides Tikhonov regularisation
  // so the solve is robust when the system is near-singular.
  Vec solve_symmetric(Mat S, Vec b, int n, double eps = 1e-9)
  {
    // Regularise
    for (int i = 0; i < n; ++i) S[i * n + i] += eps;

    // Cholesky: S = L * L^T
    Mat L = mat_zeros(n, n);
    for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j <= i; ++j)
      {
        double s = S[i * n + j];
        for (int k = 0; k < j; ++k)
          s -= L[i * n + k] * L[j * n + k];

        if (i == j)
        {
          if (s < 0.0) s = 0.0;   // clamp numerical noise
          L[i * n + i] = std::sqrt(s);
        }
        else
        {
          L[i * n + j] = (L[j * n + j] > 0.0) ? s / L[j * n + j] : 0.0;
        }
      }
    }

    // Forward substitution: L z = b
    Vec z(n, 0.0);
    for (int i = 0; i < n; ++i)
    {
      double s = b[i];
      for (int k = 0; k < i; ++k) s -= L[i * n + k] * z[k];
      z[i] = (L[i * n + i] > 0.0) ? s / L[i * n + i] : 0.0;
    }

    // Back substitution: L^T x = z
    Vec x(n, 0.0);
    for (int i = n - 1; i >= 0; --i)
    {
      double s = z[i];
      for (int k = i + 1; k < n; ++k) s -= L[k * n + i] * x[k];
      x[i] = (L[i * n + i] > 0.0) ? s / L[i * n + i] : 0.0;
    }

    return x;
  }

}

namespace arbor::utils::fitting
{
FittingComplex::FittingComplex(int K, double step_deg) : m_K(K), m_step_deg(step_deg) {}

FittingResult FittingComplex::fit(const std::vector<Vec3>& points, double tolerance)
{
  FittingResult result;
  result.shape_type = "fourier_complex";
  if (points.size() < 3) return result;

  double zsum = 0;
  for (const auto& v : points) zsum += v.z;
  m_zmean = zsum / points.size();

  result.center = calculate_centroid(points);
  PolarData polar = to_polar(points, result.center);
  double fill_radius = calculate_median(polar.r);
  PolarData augmented = inject_missing_angles(polar, fill_radius);

  Vec coefficients = fit_fourier(augmented.theta, augmented.r);

  for (size_t i = 0; i < polar.theta.size(); ++i)
  {
    if (std::abs(polar.r[i] - evaluate_fourier(polar.theta[i], coefficients)) < tolerance)
      result.inlier_indices.push_back(static_cast<int>(i));
  }

  result.success = !result.inlier_indices.empty();
  result.inlier_percentage = (100.0 * result.inlier_indices.size()) / points.size();
  result.arc_coverage_deg = calculate_arc_coverage(points, result.inlier_indices, result.center);
  if (result.arc_coverage_deg < 300) result.success = false;

  for (double c : coefficients) result.parameters.push_back(c);

  for (int i = 0; i <= 360; i += 2)
  {
    double angle_rad = i * M_PI / 180.0;
    double radius    = evaluate_fourier(angle_rad, coefficients);
    result.nodes.push_back({
      result.center.x + radius * std::cos(angle_rad),
      result.center.y + radius * std::sin(angle_rad),
      m_zmean
    });
  }

  size_t n     = result.nodes.size();
  double A     = 0.0;
  double Cx    = 0.0;
  double Cy    = 0.0;
  double ref_x = result.center.x;
  double ref_y = result.center.y;

  for (size_t i = 0; i < n; ++i)
  {
    const auto& p0 = result.nodes[i];
    const auto& p1 = result.nodes[(i + 1) % n];

    double x0 = p0.x - ref_x,  y0 = p0.y - ref_y;
    double x1 = p1.x - ref_x,  y1 = p1.y - ref_y;

    double cross = x0 * y1 - x1 * y0;
    A  += cross;
    Cx += (x0 + x1) * cross;
    Cy += (y0 + y1) * cross;
  }
  A *= 0.5;

  if (std::abs(A) < 1e-9)
  {
    result.radius = 0.0;
  }
  else
  {
    result.radius = std::sqrt(std::abs(A) / M_PI);
    Cx /= (6.0 * A);
    Cy /= (6.0 * A);
    result.center = {Cx + ref_x, Cy + ref_y, m_zmean};
  }

  return result;
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

Vec3 FittingComplex::calculate_centroid(const std::vector<Vec3>& points) const
{
  double sx = 0, sy = 0, sz = 0;
  for (const auto& p : points) { sx += p.x; sy += p.y; sz += p.z; }
  double n = static_cast<double>(points.size());
  return {sx/n, sy/n, sz/n};
}

PolarData FittingComplex::to_polar(const std::vector<Vec3>& points, const Vec3& center) const
{
  PolarData p;
  for (const auto& pt : points) {
    double dx = pt.x - center.x, dy = pt.y - center.y;
    p.theta.push_back(std::fmod(std::atan2(dy, dx) + 2.0 * M_PI, 2.0 * M_PI));
    p.r.push_back(std::sqrt(dx*dx + dy*dy));
  }
  return p;
}

double FittingComplex::calculate_median(std::vector<double> v) const
{
  if (v.empty()) return 0.0;
  std::nth_element(v.begin(), v.begin() + v.size()/2, v.end());
  return v[v.size()/2];
}

PolarData FittingComplex::inject_missing_angles(const PolarData& polar, double fill_radius) const
{
  double step_rad = m_step_deg * M_PI / 180.0;

  std::vector<size_t> idx(polar.theta.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b){
    return polar.theta[a] < polar.theta[b];
  });

  PolarData res;
  for (size_t i = 0; i < idx.size(); ++i)
  {
    size_t curr = idx[i];
    size_t next = idx[(i + 1) % idx.size()];
    res.theta.push_back(polar.theta[curr]);
    res.r.push_back(polar.r[curr]);

    bool last = (i == idx.size() - 1);
    double diff = last
    ? (2.0*M_PI - polar.theta[curr] + polar.theta[next])
      : (polar.theta[next] - polar.theta[curr]);

    if (diff > step_rad)
    {
      double end = last ? 2.0*M_PI + polar.theta[next] : polar.theta[next];
      for (double a = polar.theta[curr] + step_rad; a < end; a += step_rad) {
        res.theta.push_back(std::fmod(a, 2.0 * M_PI));
        res.r.push_back(fill_radius);
      }
    }
  }
  return res;
}

// Build the Fourier design matrix (row-major, m x (1+2K))
// col 0      : 1
// col 2k-1   : cos(k*theta)
// col 2k     : sin(k*theta)
std::vector<double> FittingComplex::build_fourier_matrix(const std::vector<double>& theta) const
{
  int m   = static_cast<int>(theta.size());
  int n   = 1 + 2 * m_K;
  Mat X   = mat_zeros(m, n);

  for (int i = 0; i < m; ++i)
  {
    X[i * n + 0] = 1.0;
    for (int k = 1; k <= m_K; ++k)
    {
      X[i * n + (2*k - 1)] = std::cos(k * theta[i]);
      X[i * n + (2*k    )] = std::sin(k * theta[i]);
    }
  }
  return X;
}

std::vector<double> FittingComplex::fit_fourier(const std::vector<double>& theta, const std::vector<double>& r) const
{
  int m = static_cast<int>(theta.size());
  int n = 1 + 2 * m_K;

  Mat X   = build_fourier_matrix(theta);
  Mat XtX = mat_ata(X, m, n);
  Vec Xty = mat_aty(X, r, m, n);

  return solve_symmetric(XtX, Xty, n);
}

double FittingComplex::evaluate_fourier(double theta, const std::vector<double>& coeffs) const
{
  double r = coeffs[0];
  for (int k = 1; k <= m_K; ++k)
    r += coeffs[2*k - 1] * std::cos(k * theta)
    + coeffs[2*k    ] * std::sin(k * theta);
  return r;
}

} // namespace arbor::utils::fitting

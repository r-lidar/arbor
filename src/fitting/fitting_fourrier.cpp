#include "fitting.h"
#include <numeric>

namespace arbor::utils::fitting {

FittingComplex::FittingComplex(int K, double step_deg) : m_K(K), m_step_deg(step_deg) {}

FittingResult FittingComplex::fit(const std::vector<Vec3>& points, double tolerance)
{
  FittingResult result;
  result.shape_type = "fourier_complex";
  if (points.size() < 3) return result;

  double zsum = 0;
  for (const auto& v : points) zsum += v.z;
  m_zmean = zsum/points.size();

  result.center = calculate_centroid(points);
  PolarData polar = to_polar(points, result.center);
  double fill_radius = calculate_median(polar.r);
  PolarData augmented = inject_missing_angles(polar, fill_radius);

  Eigen::VectorXd coefficients = fit_fourier(augmented.theta, augmented.r);

  for (size_t i = 0; i < polar.theta.size(); ++i)
  {
    if (std::abs(polar.r[i] - evaluate_fourier(polar.theta[i], coefficients)) < tolerance)
    {
      result.inlier_indices.push_back(static_cast<int>(i));
    }
  }

  result.success = !result.inlier_indices.empty();
  result.inlier_percentage = (100.0 * result.inlier_indices.size()) / points.size();
  result.arc_coverage_deg = calculate_arc_coverage(points, result.inlier_indices, result.center);
  if (result.arc_coverage_deg < 300) result.success = false;

  for (int i = 0; i < coefficients.size(); ++i) result.parameters.push_back(coefficients(i));

  for (int i = 0 ; i <= 360 ; i+=2)
  {
    double angle_deg = static_cast<double>(i);
    double angle_rad = angle_deg*M_PI/180.0;
    double radius = evaluate_fourier(angle_rad, coefficients);
    double x = result.center.x + radius*std::cos(angle_rad);
    double y = result.center.y + radius*std::sin(angle_rad);
    double z = m_zmean;
    result.nodes.push_back({x, y, z});
  }

  size_t n = result.nodes.size();

  double A = 0.0;
  double Cx = 0.0;
  double Cy = 0.0;

  for(size_t i = 0; i < n; ++i)
  {
    const auto& p0 = result.nodes[i];
    const auto& p1 = result.nodes[(i + 1) % n];

    double x0 = p0.x;
    double y0 = p0.y;
    double x1 = p1.x;
    double y1 = p1.y;

    double cross = x0 * y1 - x1 * y0;

    A  += cross;
    Cx += (x0 + x1) * cross;
    Cy += (y0 + y1) * cross;
  }

  A *= 0.5;

  Cx /= (6.0 * A);
  Cy /= (6.0 * A);

  result.center = {Cx, Cy, m_zmean};
  result.radius = std::sqrt(A/M_PI);


  return result;
}

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
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b){ return polar.theta[a] < polar.theta[b]; });

  PolarData res;
  for (size_t i = 0; i < idx.size(); ++i)
  {
    size_t curr = idx[i], next = idx[(i + 1) % idx.size()];
    res.theta.push_back(polar.theta[curr]);
    res.r.push_back(polar.r[curr]);

    double diff = (i == idx.size() - 1) ? (2.0*M_PI - polar.theta[curr] + polar.theta[next]) : (polar.theta[next] - polar.theta[curr]);
    if (diff > step_rad)
    {
      for (double a = polar.theta[curr] + step_rad; a < (i == idx.size() - 1 ? 2.0*M_PI + polar.theta[next] : polar.theta[next]); a += step_rad) {
        res.theta.push_back(std::fmod(a, 2.0*M_PI));
        res.r.push_back(fill_radius);
      }
    }
  }
  return res;
}

Eigen::MatrixXd FittingComplex::build_fourier_matrix(const std::vector<double>& theta) const
{
  Eigen::MatrixXd X(theta.size(), 1 + 2 * m_K);
  X.col(0).setOnes();
  for (int k = 1; k <= m_K; ++k)
  {
    for (size_t i = 0; i < theta.size(); ++i)
    {
      X(i, 2*k-1) = std::cos(k * theta[i]);
      X(i, 2*k) = std::sin(k * theta[i]);
    }
  }
  return X;
}

Eigen::VectorXd FittingComplex::fit_fourier(const std::vector<double>& theta, const std::vector<double>& r) const
{
  Eigen::MatrixXd X = build_fourier_matrix(theta);
  Eigen::VectorXd y = Eigen::Map<const Eigen::VectorXd>(r.data(), r.size());
  return (X.transpose() * X).ldlt().solve(X.transpose() * y);
}

double FittingComplex::evaluate_fourier(double theta, const Eigen::VectorXd& coeffs) const
{
  double r = coeffs(0);
  for (int k = 1; k <= m_K; ++k)
  {
    r += coeffs(2*k-1) * std::cos(k * theta) + coeffs(2*k) * std::sin(k * theta);
  }
  return r;
}

}

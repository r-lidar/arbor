#include "fitting.h"

namespace arbor::utils::fitting {

FittingEllipse::FittingEllipse(int max_iterations, int min_inliers, unsigned seed) : m_max_iterations(max_iterations), m_min_inliers(min_inliers), m_rng(seed)
{

}

Eigen::MatrixXd FittingEllipse::project_to_2d(const std::vector<Vec3>& points) const
{
  Eigen::MatrixXd pts(points.size(), 2);
  for (size_t i = 0; i < points.size(); ++i) {
    pts(i, 0) = points[i].x; pts(i, 1) = points[i].y;
  }
  return pts;
}

FittingEllipse::EllipseParams FittingEllipse::fit_ellipse_algebraic(const std::vector<Vec3>& pts) const
{
  EllipseParams res;
  if (pts.size() < 5) return res;
  Eigen::MatrixXd points_2d = project_to_2d(pts);
  Eigen::MatrixXd D(points_2d.rows(), 6);
  for (int i = 0; i < points_2d.rows(); ++i)
  {
    double x = points_2d(i,0), y = points_2d(i,1);
    D.row(i) << x*x, x*y, y*y, x, y, 1.0;
  }
  Eigen::VectorXd p = Eigen::JacobiSVD<Eigen::MatrixXd>(D.transpose()*D, Eigen::ComputeFullV).matrixV().col(5);
  res.a = p(0); res.b = p(1); res.c = p(2); res.d = p(3); res.e = p(4); res.f = p(5);
  res.valid = (res.b*res.b - 4*res.a*res.c < 0);
  return res;
}

std::vector<double> FittingEllipse::calculate_distances(const std::vector<Vec3>& pts, const EllipseParams& p) const
{
  std::vector<double> dists;
  for (const auto& pt : pts) {
    double f = p.a*pt.x*pt.x + p.b*pt.x*pt.y + p.c*pt.y*pt.y + p.d*pt.x + p.e*pt.y + p.f;
    double gx = 2*p.a*pt.x + p.b*pt.y + p.d, gy = p.b*pt.x + 2*p.c*pt.y + p.e;
    dists.push_back(std::abs(f) / std::sqrt(gx*gx + gy*gy + 1e-12));
  }
  return dists;
}

FittingEllipse::EllipseGeometry FittingEllipse::get_ellipse_geometry(const EllipseParams& p) const
{
  EllipseGeometry g;
  double det = p.b*p.b - 4*p.a*p.c;
  if (det >= 0) return g;
  g.cx = (2*p.c*p.d - p.b*p.e)/det; g.cy = (2*p.a*p.e - p.b*p.d)/det;
  g.angle = 0.5 * std::atan2(p.b, p.a - p.c);
  double up = 2*(p.a*p.e*p.e + p.c*p.d*p.d + p.f*p.b*p.b - p.b*p.d*p.e - 4*p.a*p.c*p.f);
  double root = std::sqrt(std::pow(p.a-p.c, 2) + p.b*p.b);
  g.major = std::sqrt(std::abs(up / (det * (p.a + p.c - root))));
  g.minor = std::sqrt(std::abs(up / (det * (p.a + p.c + root))));
  if (g.minor > g.major) std::swap(g.minor, g.major);
  g.valid = true;
  return g;
}

FittingResult FittingEllipse::fit(const std::vector<Vec3>& points, double tolerance)
{
  FittingResult best;
  best.shape_type = "ellipse";
  if (points.size() < 5) return best;

  double zsum = 0;
  for (const auto& v : points) zsum += v.z;
  m_zmean = zsum/points.size();

  std::uniform_int_distribution<size_t> dist(0, points.size() - 1);
  for (int i = 0; i < m_max_iterations; ++i)
  {
    std::vector<Vec3> sample;
    while(sample.size() < 5) sample.push_back(points[dist(m_rng)]);
    auto p = fit_ellipse_algebraic(sample);
    if (!p.valid) continue;
    auto d = calculate_distances(points, p);
    std::vector<int> inliers;
    for (int j=0; j<d.size(); ++j) if (d[j] < tolerance) inliers.push_back(j);

    if (inliers.size() > best.inlier_indices.size())
    {
      auto g = get_ellipse_geometry(p);
      best.inlier_indices = inliers;
      best.success = g.major < 3*g.minor;
    }
  }

  if (best.success)
  {
    std::vector<Vec3> inlier_pts;
    for(int idx : best.inlier_indices) inlier_pts.push_back(points[idx]);
    auto final_p = fit_ellipse_algebraic(inlier_pts);
    auto g = get_ellipse_geometry(final_p);
    best.radius = std::sqrt(g.major*g.minor);
    best.center = {g.cx, g.cy, m_zmean};
    best.parameters = {final_p.a, final_p.b, final_p.c, final_p.d, final_p.e, final_p.f, g.major, g.minor, g.angle};
    best.inlier_percentage = 100.0 * best.inlier_indices.size() / points.size();
    best.arc_coverage_deg = IFittingStrategy::calculate_arc_coverage(points, best.inlier_indices, best.center);

    double a = g.major;
    double b = g.minor;
    double theta = g.angle+M_PI/2;

    double cos_t = std::cos(theta);
    double sin_t = std::sin(theta);

    for (int i = 0; i <= 360; i += 2)
    {
      double t = i * M_PI / 180.0;

      double xr = a * std::cos(t);
      double yr = b * std::sin(t);

      double x = g.cx + xr * cos_t - yr * sin_t;
      double y = g.cy + xr * sin_t + yr * cos_t;
      double z = m_zmean;

      best.nodes.push_back({x, y, z});
    }

  }

  return best;
}

}

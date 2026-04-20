#include "fitting.h"
#include <cmath>
#include <limits>
#include <algorithm>

namespace arbor::utils::fitting {

FittingCircle::FittingCircle(int max_iterations, double early_exit_ratio, unsigned seed)
  : m_max_iterations(max_iterations),
    m_early_exit_ratio(early_exit_ratio),
    m_rng(seed)
{
}

FittingResult FittingCircle::fit(const std::vector<Vec3>& points, double tolerance)
{
  FittingResult result;
  result.shape_type = "circle";

  double zsum = 0;
  for (const auto& v : points) zsum += v.z;
  m_zmean = zsum / points.size();

  if (points.size() < 3) {
    result.success = false;
    return result;
  }

  // Fit circle using RANSAC
  CircleParams circle = fit_circle_ransac(points, tolerance);

  if (!circle.valid) {
    result.success = false;
    return result;
  }

  // Find inliers
  std::vector<int> inliers = find_inliers(points, circle, tolerance);

  if (inliers.empty()) {
    result.success = false;
    return result;
  }

  // Calculate 3D center
  Vec3 center_3d = calculate_3d_center(circle, points);

  // Calculate insider percentage (points strictly inside the circle)
  int insiders = 0;
  for (const auto& p : points) {
    double dx = p.x - circle.cx;
    double dy = p.y - circle.cy;
    double dist_to_center = std::sqrt(dx * dx + dy * dy);
    if (dist_to_center < circle.radius - tolerance) {
      ++insiders;
    }
  }

  // Calculate metrics
  result.success = true;
  result.center = center_3d;
  result.radius = circle.radius;
  result.inlier_indices = inliers;
  result.inlier_percentage = (double)inliers.size() / points.size() * 100.0;
  result.insider_percentage = (double)insiders / points.size() * 100.0;
  result.arc_coverage_deg = calculate_arc_coverage(points, inliers, center_3d);
  result.parameters = {circle.cx, circle.cy, circle.radius};

  // Generate circle nodes
  double r = circle.radius;
  for (int i = 0; i <= 360; i += 2)
  {
    double t = i * M_PI / 180.0;
    double x = center_3d.x + r * std::cos(t);
    double y = center_3d.y + r * std::sin(t);
    double z = m_zmean;

    result.nodes.push_back({x, y, z});
  }

  return result;
}

FittingCircle::CircleParams FittingCircle::fit_circle_ransac(const std::vector<Vec3>& points, double tolerance) const
{
  CircleParams best_circle;
  best_circle.valid = false;

  int n = static_cast<int>(points.size());
  if (n < 3) {
    return best_circle;
  }

  std::uniform_int_distribution<int> dist(0, n - 1);

  int max_inliers = 0;
  int early_exit_threshold = static_cast<int>(m_early_exit_ratio * n);

  for (int iter = 0; iter < m_max_iterations; ++iter)
  {
    // Pick 3 unique random points
    int idx1 = dist(m_rng);
    int idx2, idx3;
    do { idx2 = dist(m_rng); } while (idx2 == idx1);
    do { idx3 = dist(m_rng); } while (idx3 == idx1 || idx3 == idx2);

    // Fit circle on these 3 points
    CircleParams circle = fit_circle_on_3_points(
      points[idx1].x, points[idx1].y,
      points[idx2].x, points[idx2].y,
      points[idx3].x, points[idx3].y
    );

    if (!circle.valid || circle.radius <= 0.0 || std::isnan(circle.radius)) {
      continue;
    }

    // Count inliers
    int inliers = 0;
    for (int j = 0; j < n; ++j)
    {
      double dist_val = point_to_circle_distance(points[j].x, points[j].y, circle);
      if (dist_val < tolerance) {
        ++inliers;
      }
    }

    // Update best model
    if (inliers > max_inliers)
    {
      max_inliers = inliers;
      best_circle = circle;

      // Early exit if we have enough inliers
      if (max_inliers >= early_exit_threshold) {
        break;
      }
    }
  }

  if (max_inliers > 0) {
    best_circle.valid = true;
  }

  return best_circle;
}

FittingCircle::CircleParams FittingCircle::fit_circle_on_3_points(
    double x1, double y1,
    double x2, double y2,
    double x3, double y3) const
{
  CircleParams result;
  result.valid = false;

  // Calculate the coefficients for the linear system
  double A = 2.0 * (x2 - x1);
  double B = 2.0 * (y2 - y1);
  double C = x2*x2 + y2*y2 - x1*x1 - y1*y1;
  double D = 2.0 * (x3 - x1);
  double E = 2.0 * (y3 - y1);
  double G = x3*x3 + y3*y3 - x1*x1 - y1*y1;

  // Solve for cx and cy using Cramer's rule
  double denominator = A * E - B * D;

  if (std::fabs(denominator) < 1e-12) {
    return result;
  }

  double cx = (C * E - B * G) / denominator;
  double cy = (A * G - C * D) / denominator;

  // Calculate the radius
  double r = std::sqrt((x1 - cx) * (x1 - cx) + (y1 - cy) * (y1 - cy));

  result.cx = cx;
  result.cy = cy;
  result.radius = r;
  result.valid = true;

  return result;
}

double FittingCircle::point_to_circle_distance(double px, double py, const CircleParams& circle) const
{
  double dx = px - circle.cx;
  double dy = py - circle.cy;
  double dist_to_center = std::sqrt(dx * dx + dy * dy);
  return std::abs(dist_to_center - circle.radius);
}

std::vector<int> FittingCircle::find_inliers(const std::vector<Vec3>& points, const CircleParams& circle, double tolerance) const
{
  std::vector<int> inliers;

  for (size_t i = 0; i < points.size(); ++i) {
    double dist = point_to_circle_distance(points[i].x, points[i].y, circle);
    if (dist <= tolerance) {
      inliers.push_back(static_cast<int>(i));
    }
  }

  return inliers;
}

Vec3 FittingCircle::calculate_3d_center(const CircleParams& circle, const std::vector<Vec3>& points) const
{
  // Use the 2D circle center directly, with average Z
  return {circle.cx, circle.cy, m_zmean};
}

} // namespace arbor::utils::fitting

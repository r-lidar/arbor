#include "fitting.h"
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/SVD>

namespace arbor::utils::fitting {

FittingResult FittingCircle::fit(const std::vector<Vec3>& points, double tolerance)
{
  FittingResult result;
  result.shape_type = "circle";

  double zsum = 0;
  for (const auto& v : points) zsum += v.z;
  m_zmean = zsum/points.size();

  if (points.size() < 3) {
    result.success = false;
    return result;
  }

  // Project 3D points to 2D plane using PCA
  Eigen::MatrixXd points_2d = project_to_2d(points);

  // Fit circle using algebraic least squares method
  CircleParams circle = fit_circle_algebraic(points);

  if (!circle.valid) {
    result.success = false;
    return result;
  }

  // Find inliers
  std::vector<int> inliers = find_inliers(points_2d, circle, tolerance);

  if (inliers.empty()) {
    result.success = false;
    return result;
  }

  // Calculate 3D center
  Vec3 center_3d = calculate_3d_center(circle, points_2d, points);

  // Calculate metrics
  result.success = true;
  result.center = center_3d;
  result.radius = circle.radius;
  result.inlier_indices = inliers;
  result.inlier_percentage = (double)inliers.size() / points.size() * 100.0;
  result.arc_coverage_deg = calculate_arc_coverage(points, inliers, center_3d);
  result.parameters = {circle.cx, circle.cy, circle.radius};

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

FittingCircle::CircleParams FittingCircle::fit_circle_algebraic(const std::vector<Vec3>& points) const
{
  CircleParams result;
  result.valid = false;

  if (points.size() < 3) {
    return result;
  }

  // Project to 2D first
  Eigen::MatrixXd points_2d = project_to_2d(points);

  const int n = points_2d.rows();

  // Build design matrix for algebraic circle fit
  // Circle equation: x^2 + y^2 + a*x + b*y + c = 0
  Eigen::MatrixXd A(n, 3);
  Eigen::VectorXd b(n);

  for (int i = 0; i < n; ++i)
  {
    double x = points_2d(i, 0);
    double y = points_2d(i, 1);
    A(i, 0) = x;
    A(i, 1) = y;
    A(i, 2) = 1.0;
    b(i) = -(x * x + y * y);
  }

  // Solve using least squares: A * [a, b, c]^T = b
  Eigen::VectorXd solution = A.colPivHouseholderQr().solve(b);

  double a = solution(0);
  double b_coef = solution(1);
  double c = solution(2);

  // Convert to center-radius form
  result.cx = -a / 2.0;
  result.cy = -b_coef / 2.0;
  double radius_squared = (a * a + b_coef * b_coef) / 4.0 - c;

  if (radius_squared > 0) {
    result.radius = std::sqrt(radius_squared);
    result.valid = true;
  }

  return result;
}

Eigen::MatrixXd FittingCircle::project_to_2d(const std::vector<Vec3>& points) const
{
  const int n = points.size();

  // Build point matrix
  Eigen::MatrixXd P(n, 3);
  for (int i = 0; i < n; ++i) {
    P(i, 0) = points[i].x;
    P(i, 1) = points[i].y;
    P(i, 2) = points[i].z;
  }

  // Center the points
  Eigen::Vector3d centroid = P.colwise().mean();
  Eigen::MatrixXd centered = P.rowwise() - centroid.transpose();

  // Compute covariance matrix
  Eigen::Matrix3d cov = (centered.transpose() * centered) / (n - 1);

  // SVD to find principal components
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(cov, Eigen::ComputeFullU);
  Eigen::Matrix3d U = svd.matrixU();

  // Project onto first two principal components
  Eigen::MatrixXd points_2d = centered * U.leftCols(2);

  return points_2d;
}

double FittingCircle::point_to_circle_distance(double px, double py, const CircleParams& circle) const
{
  double dx = px - circle.cx;
  double dy = py - circle.cy;
  double dist_to_center = std::sqrt(dx * dx + dy * dy);
  return std::abs(dist_to_center - circle.radius);
}

std::vector<int> FittingCircle::find_inliers(const Eigen::MatrixXd& points_2d, const CircleParams& circle, double tolerance) const
{
  std::vector<int> inliers;

  for (int i = 0; i < points_2d.rows(); ++i) {
    double dist = point_to_circle_distance(points_2d(i, 0), points_2d(i, 1), circle);
    if (dist <= tolerance) {
      inliers.push_back(i);
    }
  }

  return inliers;
}

Vec3 FittingCircle::calculate_3d_center(const CircleParams& circle, const Eigen::MatrixXd& points_2d, const std::vector<Vec3>& points_3d) const
{
  const int n = points_3d.size();

  // Build point matrix and center
  Eigen::MatrixXd P(n, 3);
  for (int i = 0; i < n; ++i) {
    P(i, 0) = points_3d[i].x;
    P(i, 1) = points_3d[i].y;
    P(i, 2) = points_3d[i].z;
  }

  Eigen::Vector3d centroid = P.colwise().mean();
  Eigen::MatrixXd centered = P.rowwise() - centroid.transpose();

  // Compute covariance and get principal components
  Eigen::Matrix3d cov = (centered.transpose() * centered) / (n - 1);
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(cov, Eigen::ComputeFullU);
  Eigen::Matrix3d U = svd.matrixU();

  // Transform 2D circle center back to 3D
  Eigen::Vector2d center_2d(circle.cx, circle.cy);
  Eigen::Vector3d center_in_pca = U.leftCols(2) * center_2d;
  Eigen::Vector3d center_3d = center_in_pca + centroid;

  return {center_3d(0), center_3d(1), m_zmean};
}

} // namespace arbor::utils::fitting

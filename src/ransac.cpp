#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>

using namespace Rcpp;

// Function to fit a circle to three points using Cramer's rule
// [[Rcpp::export]]
std::vector<double> fit_circle_on_3_points_stl(const std::vector<std::vector<double>>& points_subset) {
  double x1 = points_subset[0][0], y1 = points_subset[0][1];
  double x2 = points_subset[1][0], y2 = points_subset[1][1];
  double x3 = points_subset[2][0], y3 = points_subset[2][1];

  double A = 2 * (x2 - x1), B = 2 * (y2 - y1);
  double C = x2*x2 + y2*y2 - x1*x1 - y1*y1;
  double D = 2 * (x3 - x1), E = 2 * (y3 - y1);
  double G = x3*x3 + y3*y3 - x1*x1 - y1*y1;

  double denominator = A * E - B * D;
  if (denominator == 0) return {0, 0, 0};

  double a = (C * E - B * G) / denominator;
  double b = (A * G - C * D) / denominator;
  double r = std::sqrt((x1 - a) * (x1 - a) + (y1 - b) * (y1 - b));

  return {a, b, r};
}

// Function to fit a circle to a point cloud using RANSAC
// [[Rcpp::export]]
List fit_circle_stl(NumericMatrix points, int num_iterations = 100, double inlier_threshold = 0.01) {
  int n = points.nrow();
  if (points.ncol() < 3 || n < 3) stop("Input must be an Nx3 matrix.");

  std::vector<std::vector<double>> point_vec(n, std::vector<double>(3));
  for (int i = 0; i < n; ++i) {
    point_vec[i] = {points(i, 0), points(i, 1), points(i, 2)};
  }

  std::vector<double> best_circle(3);
  int max_inliers = 0;
  std::vector<double> z_vals(n);
  std::transform(point_vec.begin(), point_vec.end(), z_vals.begin(), [](const std::vector<double>& p) { return p[2]; });

  std::random_device rd;
  std::mt19937 gen(rd());

  for (int i = 0; i < num_iterations; ++i) {
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), gen);

    std::vector<std::vector<double>> points_subset = {point_vec[indices[0]], point_vec[indices[1]], point_vec[indices[2]]};
    std::vector<double> params = fit_circle_on_3_points_stl(points_subset);
    double cx = params[0], cy = params[1], r = params[2];

    std::vector<double> distances(n);
    std::transform(point_vec.begin(), point_vec.end(), distances.begin(), [cx, cy](const std::vector<double>& p) {
      return std::sqrt((p[0] - cx) * (p[0] - cx) + (p[1] - cy) * (p[1] - cy));
    });

    std::vector<double> residuals(n);
    std::transform(distances.begin(), distances.end(), residuals.begin(), [r](double d) { return std::abs(d - r); });

    int inliers = std::count_if(residuals.begin(), residuals.end(), [inlier_threshold](double res) { return res < inlier_threshold; });

    if (inliers > max_inliers) {
      max_inliers = inliers;
      best_circle = params;
    }
  }

  double center_x = best_circle[0], center_y = best_circle[1], radius = best_circle[2];

  std::vector<double> distances(n);
  std::transform(point_vec.begin(), point_vec.end(), distances.begin(), [center_x, center_y](const std::vector<double>& p) {
    return std::sqrt((p[0] - center_x) * (p[0] - center_x) + (p[1] - center_y) * (p[1] - center_y));
  });

  std::vector<double> residuals(n);
  std::transform(distances.begin(), distances.end(), residuals.begin(), [radius](double d) { return std::abs(d - radius); });

  std::vector<int> inlier_indices;
  for (int j = 0; j < n; ++j) {
    if (residuals[j] < inlier_threshold) inlier_indices.push_back(j + 1);
  }

  double rmse = std::sqrt(std::accumulate(residuals.begin(), residuals.end(), 0.0, [](double acc, double val) {
    return acc + val * val;
  }) / n);

  // Compute angular range
  std::vector<double> angles;
  for (int j : inlier_indices) {
    angles.push_back(std::atan2(point_vec[j - 1][1] - center_y, point_vec[j - 1][0] - center_x));
  }

  for (double &angle : angles) {
    if (angle < 0) angle += 2 * M_PI;
  }

  std::sort(angles.begin(), angles.end());
  std::transform(angles.begin(), angles.end(), angles.begin(), [](double angle) { return angle * 180.0 / M_PI; });

  angles.erase(std::unique(angles.begin(), angles.end()), angles.end());
  double angle_range_degrees = 0;
  for (size_t i = 1; i < angles.size(); ++i) {
    if (angles[i] - angles[i - 1] <= 3.0) angle_range_degrees += 3.0;
  }

  return List::create(
    _["center_x"] = center_x,
    _["center_y"] = center_y,
    _["radius"] = radius,
    _["z"] = std::accumulate(z_vals.begin(), z_vals.end(), 0.0) / z_vals.size(),
    _["rmse"] = rmse,
    _["angle_range"] = angle_range_degrees,
    _["inliers"] = inlier_indices
  );
}

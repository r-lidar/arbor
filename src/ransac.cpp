#include <Rcpp.h>
#include <random>
#include <array>
#include <set>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// Define rounding function
#define ROUNDANY(x, m) round((x) / m) * m

// Optimized fit_circle_on_3_points_cpp
std::array<double, 3> fit_circle_on_3_points_cpp(double x1, double y1, double x2, double y2, double x3, double y3)
{
  double A = 2 * (x2 - x1), B = 2 * (y2 - y1);
  double C = x2 * x2 + y2 * y2 - x1 * x1 - y1 * y1;
  double D = 2 * (x3 - x1), E = 2 * (y3 - y1);
  double G = x3 * x3 + y3 * y3 - x1 * x1 - y1 * y1;

  double denominator = A * E - B * D;
  if (denominator == 0) return {0, 0, 0};

  double a = (C * E - B * G) / denominator;
  double b = (A * G - C * D) / denominator;
  double r = std::sqrt((x1 - a) * (x1 - a) + (y1 - b) * (y1 - b));

  return {a, b, r};
}

// [[Rcpp::export]]
List ransac_circle_cpp(NumericMatrix points, int num_iterations = 100, double inlier_threshold = 0.01, double early_exit = 0.8)
{
  int n = points.nrow();
  if (n < 3) stop("At least 3 points required.");

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(0, n - 1);

  std::array<double, 3> best_circle = {0, 0, 0};
  int max_inliers = 0;
  std::vector<int> inlier_indices;
  double sum_z = 0.0;

  for (int i = 0; i < n; ++i) sum_z += points(i, 2);
  double avg_z = sum_z / n;

  int early_exit_threshold = static_cast<int>(early_exit * n); // 80% inliers threshold

  std::vector<int> temp_inliers;
  temp_inliers.reserve(n);

  for (int i = 0; i < num_iterations; ++i)
  {
    temp_inliers.clear();

    // Select three distinct points
    int idx1 = dist(gen), idx2, idx3;
    do { idx2 = dist(gen); } while (idx2 == idx1);
    do { idx3 = dist(gen); } while (idx3 == idx1 || idx3 == idx2);

    // Fit a circle to these three points
    auto circle = fit_circle_on_3_points_cpp(
      points(idx1, 0), points(idx1, 1),
      points(idx2, 0), points(idx2, 1),
      points(idx3, 0), points(idx3, 1)
    );

    double cx = circle[0], cy = circle[1], r = circle[2];
    if (r == 0) continue;

    int inliers = 0;


    for (int j = 0; j < n; ++j)
    {
      double dx = points(j, 0) - cx, dy = points(j, 1) - cy;
      double dist2 = dx * dx + dy * dy;
      double res = std::abs(std::sqrt(dist2) - r);

      if (res < inlier_threshold)
      {
        temp_inliers.push_back(j + 1); // Store 1-based index
        ++inliers;
      }
    }

    if (inliers > max_inliers)
    {
      max_inliers = inliers;
      best_circle = circle;
      inlier_indices.swap(temp_inliers);

      // Early exit if 80% of points are inliers
      //if (max_inliers >= early_exit_threshold) break;
    }
  }

  double center_x = best_circle[0], center_y = best_circle[1], radius = best_circle[2];

  // Compute angular range efficiently
  std::set<int> unique_bins;
  for (int j : inlier_indices)
  {
    double angle = std::atan2(points(j - 1, 1) - center_y, points(j - 1, 0) - center_x) * 180 / M_PI;
    unique_bins.insert(static_cast<int>(ROUNDANY(angle, 3)));
  }

  return List::create(
    _["center_x"] = center_x,
    _["center_y"] = center_y,
    _["radius"] = radius,
    _["z"] = avg_z,
    _["covered_arc_degree"] = unique_bins.size() * 3,
    _["percentage_inlier"] = (double)max_inliers / n,
    _["inliers"] = inlier_indices
  );
}

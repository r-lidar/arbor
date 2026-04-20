#include "ransac.h"
#include <cmath>
#include <limits>
#include <random>
#include <set>
#include <algorithm>
#include <stdexcept>

// ---- Constructor ----
RansacCircle::RansacCircle(int iterations, double inlier_thr, double early_exit)
  : inlier_threshold(inlier_thr),
    early_exit_ratio(early_exit),
    num_iterations(iterations)
{
}

// ---- Add a point ----
void RansacCircle::add_point(double x, double y, double z)
{
  points.push_back({x, y, z});
}

// ---- Main RANSAC ----
void RansacCircle::find_circle()
{
  int n = static_cast<int>(points.size());
  if (n < 3) throw std::runtime_error("At least 3 points required to fit a circle.");

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(0, n - 1);

  Vec3 best_circle = {0, 0, 0};
  max_inliers = 0;
  max_insiders = 0;
  inlier_indices.clear();

  // Average Z
  avg_z = 0.0; for (const auto& p : points) avg_z += p.z; avg_z /= n;

  int early_exit_threshold = static_cast<int>(early_exit_ratio * n);

  std::vector<int> temp_inliers;
  temp_inliers.reserve(n);

  for (int i = 0; i < num_iterations; ++i)
  {
    temp_inliers.clear();

    // Pick 3 unique random points
    int idx1 = dist(gen), idx2, idx3;
    do { idx2 = dist(gen); } while (idx2 == idx1);
    do { idx3 = dist(gen); } while (idx3 == idx1 || idx3 == idx2);

    auto circle = fit_circle_on_3_points(
      points[idx1].x, points[idx1].y,
      points[idx2].x, points[idx2].y,
      points[idx3].x, points[idx3].y);

    double cx = circle.x, cy = circle.y, r = circle.z;
    if (r <= 0.0 || std::isnan(r)) continue;

    int inliers = 0;
    int insiders = 0;

    for (int j = 0; j < n; ++j)
    {
      double dx = points[j].x - cx;
      double dy = points[j].y - cy;
      double dist_to_center = std::sqrt(dx * dx + dy * dy);
      double res = std::abs(dist_to_center - r);

      if (res < inlier_threshold)
      {
        temp_inliers.push_back(j);
        ++inliers;
      }

      if (dist_to_center < r - inlier_threshold)
      {
        ++insiders;
      }
    }

    if (inliers > max_inliers)
    {
      max_inliers = inliers;
      max_insiders = insiders;
      best_circle = circle;
      inlier_indices = temp_inliers;

      if (max_inliers >= early_exit_threshold)
        break;
    }
  }

  center_x = best_circle.x;
  center_y = best_circle.y;
  radius   = best_circle.z;
}

// --- Find_circle with rotation ---
void RansacCircle::find_circle(const std::array<double, 3>& axis_start, const std::array<double, 3>& axis_end)
{
  // Compute rotation matrix to align normal with Z-axis
  auto R = compute_rotation_matrix(axis_start, axis_end);

  // Rotate all points
  std::vector<Vec3> rotated_points;
  rotated_points.reserve(points.size());
  for (const auto& p : points)
    rotated_points.push_back(apply_rotation(p, R));

  // Temporarily replace internal points, run RANSAC
  auto backup = points;
  points = rotated_points;
  find_circle();  // call existing method
  points = backup; // restore original points

  // Rotate center back to original 3D frame
  Vec3 center_rotated = {center_x, center_y, avg_z};

  // Compute inverse rotation (transpose of rotation matrix)
  std::array<std::array<double, 3>, 3> R_inv;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      R_inv[i][j] = R[j][i];

  Vec3 center_orig = apply_rotation(center_rotated, R_inv);

  center_x = center_orig.x;
  center_y = center_orig.y;
  avg_z    = center_orig.z;
}

// ---- Getters ----
std::array<double, 3> RansacCircle::get_center() const
{
  return {center_x, center_y, avg_z};
}

double RansacCircle::get_radius() const
{
  return radius;
}

double RansacCircle::get_inlier_percentage() const
{
  if (points.empty()) return 0.0;
  return static_cast<double>(max_inliers) / static_cast<double>(points.size());
}

double RansacCircle::get_inside_percentage() const
{
  if (points.empty()) return 0.0;
  return static_cast<double>(max_insiders) / static_cast<double>(points.size());
}

const std::vector<int>& RansacCircle::get_inliers() const
{
  return inlier_indices;
}

// ---- Arc coverage ----
double RansacCircle::get_arc_coverage() const
{
  if (inlier_indices.empty()) return 0.0;

  const double bin_size = 10.0;
  std::vector<double> angles;
  angles.reserve(inlier_indices.size());

  for (int idx : inlier_indices)
  {
    const auto& p = points[idx];
    double angle = std::atan2(p.y - center_y, p.x - center_x) * 180.0 / M_PI;
    if (angle < 0.0) angle += 360.0;
    angles.push_back(angle);
  }

  std::sort(angles.begin(), angles.end());

  std::set<int> unique_bins;
  for (double angle : angles)
  {
    int bin = static_cast<int>(std::round(angle / bin_size) * bin_size);
    if (bin == 360)  bin = 0;
    unique_bins.insert(bin);
  }
  std::vector<int> sorted_bins(unique_bins.begin(), unique_bins.end());

  if (sorted_bins.empty()) return 0.0;
  if (sorted_bins.size() == 1) return bin_size;

  int consecutive_count = 0;
  for (size_t i = 1; i < sorted_bins.size(); ++i)
  {
    if (sorted_bins[i] - sorted_bins[i-1] <= static_cast<int>(bin_size))
    {
      ++consecutive_count;
    }
  }

  // Check wraparound: is last bin close to first bin + 360?
  if (sorted_bins.back() >= 360.0 - bin_size && sorted_bins.front() <= bin_size)
  {
    ++consecutive_count;
  }

  return consecutive_count * bin_size;
}

// ---- Validation ----
bool RansacCircle::is_valid(double min_inlier_percentage, double max_inside_percentage, double min_arc_coverage_deg) const
{
  double inlier_p = get_inlier_percentage();
  double inside_p = get_inside_percentage();
  double arc = get_arc_coverage();

  if (radius <= 0.0) return false;
  if (inlier_p < min_inlier_percentage) return false;
  if (inside_p > max_inside_percentage) return false;
  if (arc < min_arc_coverage_deg) return false;

  return true;
}

// ---- Fit circle from 3 points ----
Vec3 RansacCircle::fit_circle_on_3_points(
    double x1, double y1,
    double x2, double y2,
    double x3, double y3)
{
  // Calculate the coefficients for the linear system
  double A = 2.0 * (x2 - x1);
  double B = 2.0 * (y2 - y1);
  double C = x2*x2 + y2*y2 - x1*x1 - y1*y1;
  double D = 2.0 * (x3 - x1);
  double E = 2.0 * (y3 - y1);
  double G = x3*x3 + y3*y3 - x1*x1 - y1*y1;

  // Solve for a and b using Cramer's rule
  double denominator = A * E - B * D;

  if (std::fabs(denominator) < 1e-12)
    return {0, 0, 0};

  double a = (C * E - B * G) / denominator;
  double b = (A * G - C * D) / denominator;

  // Calculate the radius
  double r = std::sqrt((x1 - a) * (x1 - a) + (y1 - b) * (y1 - b));

  return {a, b, r};
}

// --- Compute rotation matrix between two vectors ---
std::array<std::array<double, 3>, 3> RansacCircle::compute_rotation_matrix(
    const std::array<double, 3>& from,
    const std::array<double, 3>& to) const
{
  // Normalize input vectors
  auto norm = [](const std::array<double, 3>& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  };

  auto normalize = [&](const std::array<double, 3>& v) {
    double n = norm(v);
    if (n == 0) return std::array<double,3>{0,0,0};
    return std::array<double,3>{v[0]/n, v[1]/n, v[2]/n};
  };

  std::array<double, 3> a = normalize(from);
  std::array<double, 3> b = normalize(to);

  // Cross product and dot product
  std::array<double, 3> v = {a[1]*b[2] - a[2]*b[1],
                             a[2]*b[0] - a[0]*b[2],
                                               a[0]*b[1] - a[1]*b[0]};
  double s = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  double c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

  // Handle parallel or anti-parallel vectors
  if (s == 0.0) {
    if (c > 0) {
      return {{{1,0,0},{0,1,0},{0,0,1}}}; // identity
    } else {
      // 180° rotation: pick arbitrary perpendicular axis
      std::array<double, 3> axis = {1, 0, 0};
      if (std::fabs(a[0]) > 0.9) axis = {0, 1, 0};
      v = {a[1]*axis[2] - a[2]*axis[1],
           a[2]*axis[0] - a[0]*axis[2],
                                   a[0]*axis[1] - a[1]*axis[0]};
      s = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }
  }

  // Rodrigues' rotation formula
  std::array<std::array<double, 3>, 3> R;
  double vx = v[0]/s, vy = v[1]/s, vz = v[2]/s;
  double K[3][3] = {{0, -vz, vy}, {vz, 0, -vx}, {-vy, vx, 0}};
  double I[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      R[i][j] = I[i][j] + K[i][j]*s + K[i][0]*K[0][j]*(1-c)
      + K[i][1]*K[1][j]*(1-c) + K[i][2]*K[2][j]*(1-c);
    }
  }
  return R;
}

// --- Apply rotation matrix to point ---
Vec3 RansacCircle::apply_rotation(
    const Vec3& p,
    const std::array<std::array<double, 3>, 3>& R) const
{
  return {
  R[0][0]*p.x + R[0][1]*p.y + R[0][2]*p.z,
  R[1][0]*p.x + R[1][1]*p.y + R[1][2]*p.z,
  R[2][0]*p.x + R[2][1]*p.y + R[2][2]*p.z
};
}

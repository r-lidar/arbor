#include "fitting.h"
#include <cmath>
#include <algorithm>
#include <set>

namespace arbor::utils::fitting {

double IFittingStrategy::calculate_arc_coverage(const std::vector<Vec3>& points, const std::vector<int>& inliers,  const Vec3& center)
{
  if (inliers.empty()) return 0.0;
  std::vector<double> angles;
  for (int idx : inliers)
  {
    angles.push_back(std::atan2(points[idx].y - center.y, points[idx].x - center.x));
  }
  std::sort(angles.begin(), angles.end());

  double max_gap = 0.0;
  for (size_t i = 0; i < angles.size(); ++i)
  {
    double diff = (i == angles.size() - 1) ? (2.0 * M_PI - (angles[i] - angles[0])) : (angles[i+1] - angles[i]);
    max_gap = std::max(max_gap, diff);
  }
  return (2.0 * M_PI - max_gap) * (180.0 / M_PI);
}

FittingCircloid::FittingCircloid()
{
  R = {{
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}
    }};

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      R_inv[i][j] = R[j][i];
}

void FittingCircloid::add_point(double x, double y, double z)
{
  Vec3 p = {x, y, z};

  p.x -= m_origin.x;
  p.y -= m_origin.y;
  p.z -= m_origin.z;

  apply_rotation(p);
  m_points.push_back(p);
}

void FittingCircloid::clear() { m_points.clear(); }

void FittingCircloid::set_axe(const Vec3& from, const Vec3& to)
{
  // 1. Store 'from' as our local origin to translate points later
  m_origin = from;

  // 2. Calculate the directional vector of your point cloud
  Vec3 axis_dir = (to - from).normalized();

  // 3. Define the target axis (We want the point cloud aligned to Z)
  Vec3 z_axis(0.0, 0.0, 1.0);

  // 4. Compute rotation from the current direction to the Z-axis
  compute_rotation_matrix(axis_dir, z_axis);
}

FittingResult FittingCircloid::fit(double tolerance)
{
  FittingResult best_result;
  std::vector<std::unique_ptr<IFittingStrategy>> strategies;
  strategies.push_back(std::make_unique<FittingCircle>());
  strategies.push_back(std::make_unique<FittingEllipse>());
  strategies.push_back(std::make_unique<FittingComplex>());

  for (auto& strategy : strategies)
  {
    FittingResult current = strategy->fit(m_points, tolerance);
    if (current.success && current.inlier_percentage > best_result.inlier_percentage)
    {
      best_result = std::move(current);
    }
  }

  // 1. Reverse rotation for the center
  reverse_rotation(best_result.center);
  // 2. Undo the translation to put the center back in global 3D space
  best_result.center.x += m_origin.x;
  best_result.center.y += m_origin.y;
  best_result.center.z += m_origin.z;

  for (auto& v : best_result.nodes)
  {
    // Do the same for all nodes
    reverse_rotation(v);
    v.x += m_origin.x;
    v.y += m_origin.y;
    v.z += m_origin.z;
  }

  return best_result;
}

// --- Compute rotation matrix between two vectors ---
void FittingCircloid::compute_rotation_matrix(const Vec3& from, const Vec3& to)
{
  Vec3 a = from.normalized();
  Vec3 b = to.normalized();

  Vec3 v = a.cross(b);       // Rotation axis (scaled by sine of angle)
  double s = v.length();     // sine of the angle
  double c = a.dot(b);       // cosine of the angle

  // Handle Parallel or Anti-parallel edge cases
  if (s < 1e-9)
  {
    if (c > 0)
    {
      // Indentity matrix: R already initialized
      return;
    }
    else
    {
      // Case: Vectors are 180 degrees apart
      // Find an arbitrary perpendicular axis to rotate around
      Vec3 axis = (std::fabs(a.x) > 0.9) ? Vec3(0, 1, 0) : Vec3(1, 0, 0);
      v = a.cross(axis).normalized();
      s = 0.0; // sin(180) = 0
      c = -1.0;
      // For 180 deg, we proceed with the Rodrigues formula using c = -1
    }
  }
  else
  {
    // Normalize the axis for the skew-symmetric matrix
    v = v / s;
  }

  // Rodrigues' rotation formula components
  // R = I + (sin theta)K + (1 - cos theta)K^2
  double vx = v.x, vy = v.y, vz = v.z;
  double K[3][3] = {
    { 0,  -vz,  vy},
    { vz,  0,  -vx},
    {-vy,  vx,  0 }
  };

  // Compute final matrix R
  double one_minus_c = 1.0 - c;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      double identity = (i == j) ? 1.0 : 0.0;

      // Matrix multiplication for K^2 term
      double k_squared = K[i][0]*K[0][j] + K[i][1]*K[1][j] + K[i][2]*K[2][j];

      R[i][j] = identity + (K[i][j] * s) + (k_squared * one_minus_c);
    }
  }

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      R_inv[i][j] = R[j][i];
}


void FittingCircloid::apply_rotation(Vec3& p)
{
  double x, y, z;
  x = R[0][0]*p.x + R[0][1]*p.y + R[0][2]*p.z;
  y = R[1][0]*p.x + R[1][1]*p.y + R[1][2]*p.z;
  z = R[2][0]*p.x + R[2][1]*p.y + R[2][2]*p.z;
  p.x = x;
  p.y = y;
  p.z = z;
}

void FittingCircloid::reverse_rotation(Vec3& p)
{
  double x, y, z;
  x = R_inv[0][0]*p.x + R_inv[0][1]*p.y + R_inv[0][2]*p.z;
  y = R_inv[1][0]*p.x + R_inv[1][1]*p.y + R_inv[1][2]*p.z;
  z = R_inv[2][0]*p.x + R_inv[2][1]*p.y + R_inv[2][2]*p.z;
  p.x = x;
  p.y = y;
  p.z = z;
}

}

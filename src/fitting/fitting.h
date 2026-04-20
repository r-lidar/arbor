#pragma once


#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <random>
#include <Eigen/Dense>

namespace arbor::utils::fitting {

struct Vec3
{
  double x, y, z;
  Vec3() : x(0), y(0), z(0) {}
  Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
  double dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
  Vec3 cross(const Vec3& v) const { return { y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x}; }
  double length() const { return std::sqrt(x * x + y * y + z * z); }
  Vec3 normalized() const { double len = length(); if (len > 0) { return *this / len; } return {0, 0, 0}; }
  Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
  Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
  Vec3 operator*(double scalar) const { return {x * scalar, y * scalar, z * scalar}; }
  Vec3 operator/(double scalar) const { return {x / scalar, y / scalar, z / scalar}; }
  Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
  Vec3& operator-=(const Vec3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
};

struct FittingResult
{
  bool success = false;
  double radius = 0;
  double inlier_percentage = 0.0;
  double insider_percentage = 0.0;
  double arc_coverage_deg = 0.0;
  Vec3 center = {0.0, 0.0, 0.0};
  std::vector<int> inlier_indices;
  std::string shape_type;
  std::vector<Vec3> nodes;
  std::vector<double> parameters;
  bool is_valid(double min_inlier_percentage, double max_inside_percentage, double min_arc_coverage_deg) const
  {
    if (!success) return false;
    if (radius <= 0.0) return false;
    if (inlier_percentage < min_inlier_percentage) return false;
    if (insider_percentage > max_inside_percentage) return false;
    if (arc_coverage_deg < min_arc_coverage_deg) return false;
    return true;
  }
};

struct PolarData
{
  std::vector<double> theta;
  std::vector<double> r;
};

// --- Strategy Interface ---
class IFittingStrategy
{
public:
  virtual ~IFittingStrategy() = default;
  virtual FittingResult fit(const std::vector<Vec3>& points, double tolerance) = 0;

  static double calculate_arc_coverage(const std::vector<Vec3>& points, const std::vector<int>& inliers, const Vec3& center);
  double m_zmean;
};

// --- Fourier Fitting ---
class FittingComplex : public IFittingStrategy
{
public:
  FittingComplex(int K = 5, double step_deg = 10.0);
  FittingResult fit(const std::vector<Vec3>& points, double tolerance) override;

private:
  int m_K;
  double m_step_deg;

  Vec3 calculate_centroid(const std::vector<Vec3>& points) const;
  PolarData to_polar(const std::vector<Vec3>& points, const Vec3& center) const;
  double calculate_median(std::vector<double> values) const;
  PolarData inject_missing_angles(const PolarData& polar, double fill_radius) const;
  Eigen::MatrixXd build_fourier_matrix(const std::vector<double>& theta) const;
  Eigen::VectorXd fit_fourier(const std::vector<double>& theta, const std::vector<double>& r) const;
  double evaluate_fourier(double theta, const Eigen::VectorXd& coefficients) const;
};

// --- Circle Fitting ---
class FittingCircle : public IFittingStrategy
{
public:
  FittingCircle(int max_iterations = 1000, double early_exit_ratio = 0.9, unsigned seed = 64);
  FittingResult fit(const std::vector<Vec3>& points, double tolerance) override;

private:
  struct CircleParams
  {
    double cx, cy, radius;
    bool valid = false;
  };

  CircleParams fit_circle_ransac(const std::vector<Vec3>& points, double tolerance) const;
  CircleParams fit_circle_on_3_points(double x1, double y1, double x2, double y2, double x3, double y3) const;
  double point_to_circle_distance(double px, double py, const CircleParams& circle) const;
  std::vector<int> find_inliers(const std::vector<Vec3>& points, const CircleParams& circle, double tolerance) const;
  Vec3 calculate_3d_center(const CircleParams& circle, const std::vector<Vec3>& points) const;

  int m_max_iterations;
  double m_early_exit_ratio;
  mutable std::mt19937 m_rng;
};

// --- Ellipse Fitting ---
class FittingEllipse : public IFittingStrategy
{
public:
  FittingEllipse(int max_iterations = 1000, int min_inliers = 10, unsigned seed = 42);
  FittingResult fit(const std::vector<Vec3>& points, double tolerance) override;

private:
  struct EllipseParams
  {
    double a, b, c, d, e, f;
    bool valid = false;
  };

  struct EllipseGeometry
  {
    double cx, cy, major, minor, angle;
    bool valid = false;
  };

  EllipseParams fit_ellipse_algebraic(const std::vector<Vec3>& pts) const;
  std::vector<double> calculate_distances(const std::vector<Vec3>& pts, const EllipseParams& params) const;
  EllipseGeometry get_ellipse_geometry(const EllipseParams& params) const;
  Eigen::MatrixXd project_to_2d(const std::vector<Vec3>& points) const;

  int m_max_iterations;
  int m_min_inliers;
  mutable std::mt19937 m_rng;
};

// --- Master Orchestrator ---
class FittingCircloid
{
public:
  FittingCircloid();
  void set_axe(const Vec3& from, const Vec3& to);
  void add_point(double x, double y, double z);
  void clear();
  FittingResult fit(double tolerance = 0.01);

private:
  Vec3 m_origin;
  std::vector<Vec3> m_points;
  std::array<std::array<double, 3>, 3> R;
  std::array<std::array<double, 3>, 3> R_inv;
  void compute_rotation_matrix(const Vec3& from, const Vec3& to);
  void apply_rotation(Vec3& p);
  void reverse_rotation(Vec3& p);
};

} // namespace arbor::utils::fitting

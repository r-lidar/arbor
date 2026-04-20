#pragma once
#include <vector>
#include <array>

struct Vec3
{
  double x, y, z;
};

class RansacCircle
{
public:
  // ---- Constructors ----
  RansacCircle(int iterations = 100, double inlier_thr = 0.01, double early_exit = 0.8);

  // ---- Public methods ----
  void add_point(double x, double y, double z);
  void find_circle();
  void find_circle(const std::array<double, 3>& axis_start, const std::array<double, 3>& axis_end);

  std::array<double, 3> get_center() const;
  double get_radius() const;
  double get_inlier_percentage() const;
  double get_inside_percentage() const;
  double get_arc_coverage() const;
  const std::vector<int>& get_inliers() const;

  // ---- Validation ----
  bool is_valid(double min_inlier_percentage = 0.5,
                double max_inside_percentage = 0.3,
                double min_arc_coverage_deg = 120.0) const;

private:
  // ---- Data ----
  std::vector<Vec3> points;

  double center_x = 0.0;
  double center_y = 0.0;
  double radius = 0.0;
  double avg_z = 0.0;

  double inlier_threshold = 0.01;
  double early_exit_ratio = 0.8;
  int num_iterations = 100;

  int max_inliers = 0;
  int max_insiders = 0;
  std::vector<int> inlier_indices;

  static Vec3 fit_circle_on_3_points(double x1, double y1, double x2, double y2, double x3, double y3);
  std::array<std::array<double, 3>, 3> compute_rotation_matrix(const std::array<double, 3>& from, const std::array<double, 3>& to) const;
  Vec3 apply_rotation(const Vec3& p, const std::array<std::array<double, 3>, 3>& R) const;
};

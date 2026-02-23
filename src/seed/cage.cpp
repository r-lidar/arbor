#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <set>

using namespace Rcpp;

// Structure to hold a 3D point
struct Point3D
{
  double X;
  double Y;
  double Z;

  Point3D(double x, double y, double z) : X(x), Y(y), Z(z) {}
};

// Structure to hold a circle
struct Circle {
  double X;
  double Y;
  double Z;
  double R;
  int id;

  Circle(double x, double y, double z, double r, int i) : X(x), Y(y), Z(z), R(r), id(i) {}
};

// Generate points on a circle circumference
std::vector<Point3D> generate_circle_points(double x, double y, double z, double r, double step)
{
  std::vector<Point3D> points;

  // Calculate circumference and number of points
  double circumference = 2.0 * M_PI * r;
  int n_points = static_cast<int>(std::ceil(circumference / step));

  // Safety check
  if (n_points <= 0)
  {
    return points;
  }

  points.reserve(n_points);

  // Generate points at regular angular intervals
  for (int i = 0; i < n_points; ++i)
  {
    double theta = (2.0 * M_PI * i) / n_points;
    double px = x + r * std::cos(theta);
    double py = y + r * std::sin(theta);
    points.push_back(Point3D(px, py, z));
  }

  return points;
}

// Generate points on radii of a disk (8 points at equal angles)
std::vector<Point3D> generate_disk_radii(double x, double y, double z, double r, int n = 8)
{
  std::vector<Point3D> points;
  points.reserve(n);

  for (int i = 0; i < n; ++i) {
    double angle = (2.0 * M_PI * i) / n;
    double px = x + r * std::cos(angle);
    double py = y + r * std::sin(angle);
    points.push_back(Point3D(px, py, z));
  }

  return points;
}

// Generate connectors for a single group of circles
std::vector<Point3D> generate_connectors(std::vector<Circle> group, double step_z)
{
  std::vector<Point3D> connectors;

  if (group.size() < 2)
  {
    return connectors;
  }

  // Sort by Z coordinate
  std::sort(group.begin(), group.end(), [](const Circle& a, const Circle& b) { return a.Z < b.Z; });

  // Iterate through pairs of consecutive circles
  for (size_t i = 0; i < group.size() - 1; ++i)
  {
    const Circle& c1 = group[i];
    const Circle& c2 = group[i + 1];

    double dz = c2.Z - c1.Z;

    // Safety check for zero or negative height difference
    if (dz <= 0.0)
    {
      continue;
    }

    // Generate vertical sequence
    std::vector<double> z_seq;
    double z_current = c1.Z;

    while (z_current <= c2.Z - step_z)
    {
      z_seq.push_back(z_current);
      z_current += step_z;
    }

    // Always include the top
    z_seq.push_back(c2.Z);

    // Generate points at each z level
    for (size_t j = 0; j < z_seq.size(); ++j)
    {
      double z = z_seq[j];

      // Interpolation parameter
      double t = (z - c1.Z) / dz;

      // Interpolate center and radius
      double x_interp = c1.X + t * (c2.X - c1.X);
      double y_interp = c1.Y + t * (c2.Y - c1.Y);
      double r_interp = c1.R + t * (c2.R - c1.R);

      // Generate disk radii points
      std::vector<Point3D> disk_points = generate_disk_radii(x_interp, y_interp, z, r_interp, 8);

      // Reserve space before inserting
      connectors.reserve(connectors.size() + disk_points.size());
      connectors.insert(connectors.end(), disk_points.begin(), disk_points.end());
    }
  }

  return connectors;
}

// Generate all connectors for all groups
std::vector<Point3D> generate_all_connectors(const std::vector<Circle>& circles, double step_z)
{
  std::vector<Point3D> all_connectors;

  if (circles.empty()) {
    return all_connectors;
  }

  // Find unique IDs using set for efficiency
  std::set<int> unique_ids_set;
  for (const auto& circle : circles) {
    unique_ids_set.insert(circle.id);
  }

  std::vector<int> unique_ids(unique_ids_set.begin(), unique_ids_set.end());

  // Process each group
  for (int id : unique_ids)
  {
    // Extract circles for this group
    std::vector<Circle> group;
    for (const auto& circle : circles)
    {
      if (circle.id == id) {
        group.push_back(circle);
      }
    }

    // Generate connectors for this group
    if (group.size() > 1)
    {
      std::vector<Point3D> group_connectors = generate_connectors(group, step_z);

      // Reserve space before inserting
      all_connectors.reserve(all_connectors.size() + group_connectors.size());
      all_connectors.insert(all_connectors.end(), group_connectors.begin(), group_connectors.end());
    }
  }

  return all_connectors;
}

// Pure C++ function that does all the work
std::vector<Point3D> generate_cage(const std::vector<Circle>& circles, double decimation)
{
  // Calculate step size
  double res = decimation * 0.75;

  // Safety check
  if (res <= 0.0)
  {
    Rcpp::stop("Invalid decimation value, resulting in non-positive step size");
  }

  // Generate circle points
  std::vector<Point3D> all_circle_points;
  for (const auto& circle : circles)
  {
    std::vector<Point3D> pts = generate_circle_points(circle.X, circle.Y, circle.Z, circle.R, res);
    all_circle_points.reserve(all_circle_points.size() + pts.size());
    all_circle_points.insert(all_circle_points.end(), pts.begin(), pts.end());
  }

  // Generate connectors
  std::vector<Point3D> all_connectors = generate_all_connectors(circles, res);

  // Combine all points
  std::vector<Point3D> cage_points;
  cage_points.reserve(all_circle_points.size() + all_connectors.size());
  cage_points.insert(cage_points.end(), all_circle_points.begin(), all_circle_points.end());
  cage_points.insert(cage_points.end(), all_connectors.begin(), all_connectors.end());

  return cage_points;
}

// Rcpp wrapper - only converts input/output
// [[Rcpp::export]]
DataFrame generate_cage_cpp(DataFrame circles, double decimation)
{
  // Convert input: R DataFrame -> C++ vector of Circles
  NumericVector X = circles["X"];
  NumericVector Y = circles["Y"];
  NumericVector Z = circles["Z"];
  NumericVector R = circles["R"];
  IntegerVector id = circles["id"];

  int n = X.size();

  if (n == 0) {
    return DataFrame::create(
      Named("X") = NumericVector(),
      Named("Y") = NumericVector(),
      Named("Z") = NumericVector()
    );
  }

  std::vector<Circle> circle_vec;
  circle_vec.reserve(n);
  for (int i = 0; i < n; ++i) {

    circle_vec.push_back(Circle(X[i], Y[i], Z[i], R[i], id[i]));
  }

  // Call pure C++ function
  std::vector<Point3D> cage_points = generate_cage(circle_vec, decimation);

  // Convert output: C++ vector of Points -> R DataFrame
  size_t total_points = cage_points.size();
  NumericVector out_X(total_points);
  NumericVector out_Y(total_points);
  NumericVector out_Z(total_points);

  for (size_t i = 0; i < total_points; ++i) {
    out_X[i] = cage_points[i].X;
    out_Y[i] = cage_points[i].Y;
    out_Z[i] = cage_points[i].Z;
  }

  return DataFrame::create(
    Named("X") = out_X,
    Named("Y") = out_Y,
    Named("Z") = out_Z
  );
}

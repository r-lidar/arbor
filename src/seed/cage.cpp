#include <vector>
#include <cmath>
#include <algorithm>

struct CircleResult
{
  double X;
  double Y;
  double Z;
  double R;
  int id;

  CircleResult(double x, double y, double z, double r, int cluster_id) : X(x), Y(y), Z(z), R(r), id(cluster_id) {}
};

struct Point3D
{
  double X;
  double Y;
  double Z;

  Point3D(double x, double y, double z) : X(x), Y(y), Z(z) {}
};

struct CageResult
{
  std::vector<Point3D> points;
  // Optional: add passage and hag attributes if needed
  // std::vector<double> passage;
  // std::vector<double> hag;
};

// Generate points along a circle perimeter
std::vector<Point3D> generate_circle_points(double x, double y, double z, double r, double step = 0.1)
{
  std::vector<Point3D> points;

  // Calculate circumference and number of points
  double circumference = 2.0 * M_PI * r;
  int n_points = static_cast<int>(std::ceil(circumference / step));

  // Generate points around the circle
  for (int i = 0; i < n_points; ++i)
  {
    double theta = (2.0 * M_PI * i) / n_points;
    double px = x + r * std::cos(theta);
    double py = y + r * std::sin(theta);
    points.emplace_back(px, py, z);
  }

  return points;
}

// Generate points arranged as radii on a disk (8 directions)
std::vector<Point3D> generate_disk_radii(double x, double y, double z, double r, int n = 8)
{
  std::vector<Point3D> points;

  for (int i = 0; i < n; ++i)
  {
    double angle = (2.0 * M_PI * i) / n;
    double px = x + r * std::cos(angle);
    double py = y + r * std::sin(angle);
    points.emplace_back(px, py, z);
  }

  return points;
}

// Generate connecting points between circles in a group
struct GroupResult
{
  std::vector<Point3D> disks;
  std::vector<Point3D> centerline;
};

GroupResult generate_group_points(std::vector<CircleResult>& group, double step_z = 0.05)
{
  // Sort by Z coordinate
  std::sort(group.begin(), group.end(), [](const CircleResult& a, const CircleResult& b) {
    return a.Z < b.Z;
  });

  GroupResult result;

  for (size_t i = 0; i < group.size() - 1; ++i)
  {
    const CircleResult& c1 = group[i];
    const CircleResult& c2 = group[i + 1];

    // Generate Z sequence from c1.Z to c2.Z with step_z increments
    std::vector<double> z_seq;
    for (double z = c1.Z; z < c2.Z; z += step_z)
    {
      z_seq.push_back(z);
    }
    // Ensure the top circle is included
    if (z_seq.empty() || z_seq.back() != c2.Z)
    {
      z_seq.push_back(c2.Z);
    }

    // Interpolate for each Z level
    for (double z : z_seq)
    {
      double t = (z - c1.Z) / (c2.Z - c1.Z);

      // Interpolate center position and radius
      double x_interp = c1.X + t * (c2.X - c1.X);
      double y_interp = c1.Y + t * (c2.Y - c1.Y);
      double r_interp = c1.R + t * (c2.R - c1.R);

      // Add centerline point
      result.centerline.emplace_back(x_interp, y_interp, z);

      // Add disk radii points
      auto disk_points = generate_disk_radii(x_interp, y_interp, z, r_interp);
      result.disks.insert(result.disks.end(), disk_points.begin(), disk_points.end());
    }
  }

  return result;
}

// Generate points for all groups
GroupResult generate_all_groups(std::vector<CircleResult>& circles, double step_z = 0.05)
{
  GroupResult combined_result;

  // Find unique group IDs
  std::vector<int> unique_ids;
  for (const auto& circle : circles)
  {
    if (std::find(unique_ids.begin(), unique_ids.end(), circle.id) == unique_ids.end()) {
      unique_ids.push_back(circle.id);
    }
  }

  // Process each group
  for (int group_id : unique_ids)
  {
    // Extract circles belonging to this group
    std::vector<CircleResult> group;
    for (const auto& circle : circles)
    {
      if (circle.id == group_id)
      {
        group.push_back(circle);
      }
    }

    // Generate points for this group
    GroupResult group_result = generate_group_points(group, step_z);

    // Combine results
    combined_result.disks.insert(combined_result.disks.end(), group_result.disks.begin(), group_result.disks.end());
    combined_result.centerline.insert(combined_result.centerline.end(), group_result.centerline.begin(), group_result.centerline.end());
  }

  return combined_result;
}

// Main function: generate cage from circles
CageResult generate_cage(std::vector<CircleResult>& circles, double decimation = 0.1)
{
  CageResult cage;

  // Step 1: Generate circle perimeter points
  double res = decimation * 0.75;
  for (const auto& circle : circles)
  {
    auto circle_points = generate_circle_points(circle.X, circle.Y, circle.Z, circle.R, res);
    cage.points.insert(cage.points.end(), circle_points.begin(), circle_points.end());
  }

  // Step 2: Generate connector points between circles
  GroupResult connectors = generate_all_groups(circles, res);

  // Step 3: Combine all points
  cage.points.insert(cage.points.end(), connectors.disks.begin(), connectors.disks.end());

  return cage;
}

/**
 * @file SeedDetector_cage.cpp
 * Project: Arbor
 *
 * Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <vector>
#include <algorithm>
#include <set>

#include "SeedDetector.h"

namespace arbor::seeds {

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

std::vector<Point3D> SeedDetector::generate_cage(const std::vector<Circle>& circles, double decimation)
{
  double res = decimation * 0.75;

  if (res <= 0.0)
    throw std::runtime_error("Invalid decimation value, resulting in non-positive step size");

  std::vector<Point3D> all_circle_points;
  for (const auto& circle : circles)
  {
    std::vector<Point3D> pts = generate_circle_points(circle.X, circle.Y, circle.Z, circle.R, res);
    all_circle_points.insert(all_circle_points.end(), std::make_move_iterator(pts.begin()), std::make_move_iterator(pts.end()));
  }

  std::vector<Point3D> all_connectors = generate_all_connectors(circles, res);

  std::vector<Point3D> cage_points;
  cage_points.reserve(all_circle_points.size() + all_connectors.size());
  cage_points.insert(cage_points.end(), all_circle_points.begin(), all_circle_points.end());
  cage_points.insert(cage_points.end(), all_connectors.begin(), all_connectors.end());

  return cage_points;
}

void SeedDetector::make_cages()
{
  std::vector<Point3D> pcage = generate_cage(circles, params.pathfinder.decimation); // static function exported to R

  cages = PointCloud(pcage.size(), true);
  for (size_t i = 0 ; i < cages.size() ; i++)
  {
    cages.set_x(i, pcage[i].X);
    cages.set_y(i, pcage[i].Y);
    cages.set_z(i, pcage[i].Z);
    cages.set_passage(i, 9999); // Need to be assigned any value > 0 to be retained in find_primary_seeds()
  }
}

}

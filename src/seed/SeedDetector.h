/**
 * @file SeedDetector.h
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

#ifndef SEEDDETECTOR_H
#define SEEDDETECTOR_H

#include "arbor.h"

namespace arbor::seeds {

struct Circle
{
  double X;
  double Y;
  double Z;
  double R;
  int id;

  Circle(double x, double y, double z, double r, int i) : X(x), Y(y), Z(z), R(r), id(i) {}
};

struct Point3D
{
  double X;
  double Y;
  double Z;

  Point3D(double x, double y, double z) : X(x), Y(y), Z(z) {}
};

class SeedDetector
{
public:
  SeedDetector(const settings::ArborParameters& par) : params(par) {};
  void run(const PointCloud& scene);
  static std::vector<Circle> detect_tree_circles(const PointCloud& wood, double resolution = 0.05, int connectivity = 26, int num_ransac_iterations = 400, double inlier_threshold = 0.02, size_t min_cluster_size = 10);
  static std::vector<Point3D> generate_cage(const std::vector<Circle>& circles, double decimation);

  const PointCloud& get_long_passages() const { return long_passages; }
  const PointCloud& get_short_passages() const { return short_passages; }
  const PointCloud& get_wood() const { return wood; }
  const PointCloud& get_cages() const { return cages; }
  const PointCloud& get_primary_seeds() const { return primary_seeds; }
  const PointCloud& get_secondary_seeds() const { return secondary_seeds; }
  const PointCloud& get_seeds() const { return seeds; }
  PointCloud move_seeds() { return std::move(seeds); }   // allow to move the data out (transfer ownership)

private:
  void compute_min_hag(const PointCloud& scene);
  void slice_wood(const PointCloud& scene);
  void extract_passages(const PointCloud& scene);
  void make_cages();
  void safe_zone();
  void find_primary_seeds();
  void find_secondary_seeds();
  void merge_short_passages();
  void filter_seeds();
  PointCloud densify_passages(const PointCloud& x);


private:
  std::vector<Circle> circles;

  PointCloud long_passages;
  PointCloud short_passages;
  PointCloud wood;
  PointCloud cages;

  PointCloud primary_seeds;
  PointCloud secondary_seeds;
  PointCloud seeds;

  double min_hag;

  settings::ArborParameters params;
};

}

#endif

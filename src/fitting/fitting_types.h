/**
 * @file fitting_types.h
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

#pragma once

#include <cmath>
#include <string>
#include <vector>

namespace arbor::utils::fitting {

struct Vec3
{
  double x, y, z;
  Vec3() : x(0), y(0), z(0) {}
  Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
  double dot(const Vec3& v)   const { return x*v.x + y*v.y + z*v.z; }
  Vec3 cross(const Vec3& v)   const { return {y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x}; }
  double length()             const { return std::sqrt(x*x + y*y + z*z); }
  Vec3 normalized()           const { double len = length(); return len > 0 ? *this / len : Vec3{0,0,0}; }
  Vec3  operator+ (const Vec3& v) const { return {x+v.x, y+v.y, z+v.z}; }
  Vec3  operator- (const Vec3& v) const { return {x-v.x, y-v.y, z-v.z}; }
  Vec3  operator* (double s)      const { return {x*s,   y*s,   z*s  }; }
  Vec3  operator/ (double s)      const { return {x/s,   y/s,   z/s  }; }
  Vec3& operator+=(const Vec3& v) { x+=v.x; y+=v.y; z+=v.z; return *this; }
  Vec3& operator-=(const Vec3& v) { x-=v.x; y-=v.y; z-=v.z; return *this; }
};

struct FittingResult
{
  bool  success              = false;
  float radius               = 0.0f;
  float inlier_percentage    = 0.0f;
  float interior_percentage  = 0.0f;  // points strictly inside the fitted shape
  float arc_coverage_deg     = 0.0f;
  Vec3  center               = {0.0, 0.0, 0.0};
  std::vector<int>    inlier_indices;
  std::string         shape_type;
  std::vector<Vec3>   contour;        // outline sample points of the fitted shape
  std::vector<double> parameters;

  bool is_valid(double min_inlier_percentage,
                double max_interior_percentage,
                double min_arc_coverage_deg) const
  {
    if (!success)                                         return false;
    if (radius              <= 0.0)                       return false;
    if (inlier_percentage   < min_inlier_percentage)      return false;
    if (interior_percentage > max_interior_percentage)    return false;
    if (arc_coverage_deg    < min_arc_coverage_deg)       return false;
    return true;
  }
};

struct PolarData
{
  std::vector<double> theta;
  std::vector<double> r;
};

} // namespace arbor::utils::fitting

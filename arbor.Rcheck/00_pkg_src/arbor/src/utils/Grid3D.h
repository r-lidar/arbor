/**
 * @file Grid3D.h
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

#ifndef GRID3D_H
#define GRID3D_H

#include <vector>
#include <unordered_map>
#include <cstdint>

#include "PointCloud.h"

struct Offset
{
  int dx, dy, dz;
};

class Grid3D
{
public:
  Grid3D(const PointCloud& pc, double res);

  // Grid dimensions
  int64_t ncols, nrows, nlayers;

  // Bounding box
  double xmin, ymin, zmin;
  double xmax, ymax, zmax;

  // Inverse resolution
  double inv_xres, inv_yres, inv_zres;

  int npoints;

  std::vector<int> connected_components(int connectivity);

private:
  const PointCloud& pc;

  int64_t get_grid_index(int64_t r, int64_t c, int64_t l) const;
  void decode_grid_index(int64_t idx, int64_t& r, int64_t& c, int64_t& l) const;
  std::vector<Offset> get_neighbor_offsets(int connectivity) const;
};

#endif // GRID3D_H

/**
 * @file Raster.cpp
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

#include "Raster.h"

#include <algorithm>
#include <stdexcept>

namespace arbor {

Raster::Raster(double x_min, double y_min, double x_max, double y_max, double resolution) : res(resolution)
{
  if (!(res > 0)) throw std::runtime_error("Raster: the resolution must be positive.");

  xmin = std::floor(x_min);
  ymin = std::floor(y_min);

  // Number of cells needed to fully cover the domain
  double cols = std::ceil((x_max - xmin) / res);
  double rows = std::ceil((y_max - ymin) / res);
  if (!(cols > 0) || !(rows > 0)) throw std::runtime_error("Raster: invalid dimensions.");

  ncols = static_cast<size_t>(cols);
  nrows = static_cast<size_t>(rows);
  ymax  = ymin + nrows * res;

  values.assign(ncols * nrows, NA);
}

bool Raster::contains(double x, double y) const
{
  return !empty() && x >= xmin && x <= get_xmax() && y >= ymin && y <= ymax;
}

size_t Raster::cell_from_xy(double x, double y) const
{
  // The east and south edges belong to the last column and row
  size_t col = std::min(static_cast<size_t>((x - xmin) / res), ncols - 1);
  size_t row = std::min(static_cast<size_t>((ymax - y) / res), nrows - 1);
  return cell_from_row_col(row, col);
}

float Raster::get_value(double x, double y) const
{
  if (!contains(x, y)) return NA;
  return values[cell_from_xy(x, y)];
}

PointCloud Raster::to_points() const
{
  size_t n = values.size() - std::count_if(values.begin(), values.end(), [](float v) { return std::isnan(v); });

  PointCloud points(n);
  size_t j = 0;
  for (size_t cell = 0 ; cell < values.size() ; cell++)
  {
    if (is_na(cell)) continue;
    points.set_x(j, x_from_cell(cell));
    points.set_y(j, y_from_cell(cell));
    points.set_z(j, values[cell]);
    j++;
  }

  return points;
}

}

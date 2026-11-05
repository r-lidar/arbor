/**
 * @file Raster.h
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

#ifndef ARBOR_RASTER_H
#define ARBOR_RASTER_H

#include "PointCloud.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace arbor {

// A single band raster of square cells (e.g. a DTM). Only the values are stored (float, NaN = no data): the
// coordinates of a cell are derived from the extent, in double precision.
//
// Layout: row-major, row 0 is the north edge (ymax), column 0 the west edge (xmin). The cell coordinates are the
// centers of the cells. Same geometry as IncrementalDelaunay::Grid: xmin and ymin are floored, the extent is rounded
// up to cover the requested bounds with whole cells.
class Raster
{
public:
  static constexpr float NA = std::numeric_limits<float>::quiet_NaN();

  Raster() = default;
  Raster(double xmin, double ymin, double xmax, double ymax, double res); // every cell NA

  // Geometry
  bool   empty()      const { return values.empty(); }
  size_t get_ncols()  const { return ncols; }
  size_t get_nrows()  const { return nrows; }
  size_t get_ncells() const { return values.size(); }
  double get_res()    const { return res; }
  double get_xmin()   const { return xmin; }
  double get_ymin()   const { return ymin; }
  double get_xmax()   const { return xmin + ncols * res; }
  double get_ymax()   const { return ymax; }

  // Cells
  size_t cell_from_row_col(size_t row, size_t col) const { return row * ncols + col; }
  size_t row_from_cell(size_t cell)                const { return cell / ncols; }
  size_t col_from_cell(size_t cell)                const { return cell % ncols; }
  double x_from_col(size_t col)                    const { return xmin + (col + 0.5) * res; }
  double y_from_row(size_t row)                    const { return ymax - (row + 0.5) * res; }
  double x_from_cell(size_t cell)                  const { return x_from_col(col_from_cell(cell)); }
  double y_from_cell(size_t cell)                  const { return y_from_row(row_from_cell(cell)); }
  bool   contains(double x, double y)              const;
  size_t cell_from_xy(double x, double y)          const; // contains(x, y) must be true

  // Values
  float get_value(size_t cell)              const { return values[cell]; }
  float get_value(double x, double y)       const; // value of the cell that contains (x, y), NA outside
  void  set_value(size_t cell, float value)       { values[cell] = value; }
  bool  is_na(size_t cell)                  const { return std::isnan(values[cell]); }
  const std::vector<float> &get_values()    const { return values; }

  // One point per cell with a value (NA cells skipped), at the center of the cell, in the order of the cells.
  PointCloud to_points() const;

private:
  size_t ncols = 0;
  size_t nrows = 0;
  double xmin  = 0;
  double ymin  = 0;
  double ymax  = 0;
  double res   = 0;
  std::vector<float> values;
};

}

#endif

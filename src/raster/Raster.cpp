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

#include "myomp.h"

#include <algorithm>
#include <stdexcept>

namespace arbor {

Raster::Raster(double x_min, double y_min, double x_max, double y_max, double resolution, bool floor_origin) : res(resolution)
{
  if (!(res > 0)) throw std::runtime_error("Raster: the resolution must be positive.");

  const bool exact = !floor_origin;
  xmin = exact ? x_min : std::floor(x_min);
  ymin = exact ? y_min : std::floor(y_min);

  // Number of cells needed to fully cover the domain. An exact grid is a whole number of cells: the rounding errors of
  // its bounds must not add one.
  const double tolerance = exact ? 1e-6 : 0;
  double cols = std::ceil((x_max - xmin) / res - tolerance);
  double rows = std::ceil((y_max - ymin) / res - tolerance);
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

float Raster::nearest_value(long row, long col) const
{
  const long nr = static_cast<long>(nrows);
  const long nc = static_cast<long>(ncols);
  row = std::clamp(row, 0L, nr - 1);
  col = std::clamp(col, 0L, nc - 1);

  // Square rings of growing radius: the first ring holding values gives the nearest of them (distance in cells).
  const long max_radius = std::max(nr, nc);
  for (long radius = 0 ; radius <= max_radius ; ++radius)
  {
    float  best      = NA;
    double best_dist = std::numeric_limits<double>::max();
    for (long r = row - radius ; r <= row + radius ; ++r)
    {
      if (r < 0 || r >= nr) continue;
      const bool edge_row = r == row - radius || r == row + radius;
      const long step     = edge_row ? 1 : 2 * radius; // inside rows: only the two ends belong to the ring
      for (long c = col - radius ; c <= col + radius ; c += std::max(step, 1L))
      {
        if (c < 0 || c >= nc) continue;
        const size_t cell = cell_from_row_col(r, c);
        if (is_na(cell)) continue;
        const double dist = double(r - row) * (r - row) + double(c - col) * (c - col);
        if (dist < best_dist) { best_dist = dist; best = values[cell]; }
      }
    }
    if (!std::isnan(best)) return best;
  }
  return NA;
}

float Raster::interpolate(double x, double y) const
{
  if (empty()) return NA;

  // Position in cell centers: (0, 0) is the center of the north west cell.
  const double fc = (x - xmin) / res - 0.5;
  const double fr = (ymax - y) / res - 0.5;
  const long   nr = static_cast<long>(nrows);
  const long   nc = static_cast<long>(ncols);

  const long   c0 = static_cast<long>(std::floor(fc));
  const long   r0 = static_cast<long>(std::floor(fr));
  const double tx = std::clamp(fc - c0, 0.0, 1.0);
  const double ty = std::clamp(fr - r0, 0.0, 1.0);

  double sum = 0;
  double weights = 0;
  for (int dr = 0 ; dr <= 1 ; ++dr)
  {
    for (int dc = 0 ; dc <= 1 ; ++dc)
    {
      const long r = std::clamp(r0 + dr, 0L, nr - 1);
      const long c = std::clamp(c0 + dc, 0L, nc - 1);
      const float v = values[cell_from_row_col(r, c)];
      if (std::isnan(v)) continue;
      const double w = (dc ? tx : 1 - tx) * (dr ? ty : 1 - ty);
      sum     += w * v;
      weights += w;
    }
  }

  if (weights > 0) return static_cast<float>(sum / weights);
  return nearest_value(static_cast<long>(std::lround(fr)), static_cast<long>(std::lround(fc)));
}

void Raster::compute_hag(PointCloud &scene) const
{
  if (!scene.has_hag()) throw std::runtime_error("Raster::compute_hag(): the HAG attribute is not allocated.");
  if (std::all_of(values.begin(), values.end(), [](float v) { return std::isnan(v); }))
    throw std::runtime_error("Raster::compute_hag(): the raster has no value.");

  const long n = static_cast<long>(scene.size());

  #pragma omp parallel for
  for (long i = 0 ; i < n ; ++i)
  {
    const double ground = interpolate(scene.get_x(i), scene.get_y(i));
    scene.set_hag(i, scene.get_z(i) - ground);
  }
}

Raster Raster::focal(int window, const std::function<float(std::vector<float>&)> &fun) const
{
  Raster out = *this;
  std::fill(out.values.begin(), out.values.end(), NA);

  const long half_window = std::max(0, window / 2);
  const long nr = static_cast<long>(nrows);
  const long nc = static_cast<long>(ncols);

  #pragma omp parallel
  {
    std::vector<float> neighbors; // one buffer per thread

    #pragma omp for
    for (long row = 0 ; row < nr ; ++row)
    {
      for (long col = 0 ; col < nc ; ++col)
      {
        const size_t cell = cell_from_row_col(row, col);
        if (is_na(cell)) continue;

        const long first_row = std::max(0L, row - half_window);
        const long last_row  = std::min(nr - 1, row + half_window);
        const long first_col = std::max(0L, col - half_window);
        const long last_col  = std::min(nc - 1, col + half_window);

        neighbors.clear();
        for (long r = first_row ; r <= last_row ; ++r)
        {
          for (long c = first_col ; c <= last_col ; ++c)
          {
            const size_t neighbor = cell_from_row_col(r, c);
            if (is_na(neighbor)) continue;
            neighbors.push_back(values[neighbor]);
          }
        }

        out.values[cell] = fun(neighbors);
      }
    }
  }

  return out;
}

}

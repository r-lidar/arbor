/**
 * @file dtm.cpp
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

#include "arbor.h"
#include "hporro/delaunay.h"

namespace arbor::dtm {

PointCloud dtm(const PointCloud& scene)
{
  ServiceLocator::logger()("Computing DTM");

  size_t n = scene.size();

  // Initialize bounds
  double xmin = std::numeric_limits<double>::max(), xmax = -std::numeric_limits<double>::max();
  double ymin = std::numeric_limits<double>::max(), ymax = -std::numeric_limits<double>::max();
  double zmin = std::numeric_limits<double>::max(), zmax = -std::numeric_limits<double>::max();
  for (unsigned int i = 0; i < n ; ++i)
  {
    double x = scene.get_x(i);
    double y = scene.get_y(i);
    double z = scene.get_z(i);

    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
  }

  double x_offset = (xmin + xmax) * 0.5;
  double y_offset = (ymin + ymax) * 0.5;
  //double z_offset = (zmin + zmax) * 0.5;

  IncrementalDelaunay::Grid index(xmin - x_offset, ymin - y_offset, xmax - x_offset, ymax - y_offset, 0.5);
  IncrementalDelaunay::Triangulation d(index);


  for (unsigned int i = 0 ; i < n ; i++)
  {
    if (scene.is_ground(i))
    {
       double x = scene.get_x(i) - x_offset;
       double y = scene.get_y(i) - y_offset;
       double z = scene.get_z(i);
       IncrementalDelaunay::Vec2 p(x,y,z);
       int t = d.findContainerTriangleFast(p);
       d.delaunayInsertion(p, t);
    }
  }

  IncrementalDelaunay::Grid dtm_grid(xmin, ymin, xmax, ymax, 0.1);
  PointCloud dtm(dtm_grid.get_ncells());
  for (int i = 0 ; i < dtm_grid.get_ncells() ; i++)
  {
    double x = dtm_grid.x_from_cell(i);
    double y = dtm_grid.y_from_cell(i);
    double z = d.get_z(x - x_offset, y - y_offset);
    dtm.set_x(i, x);
    dtm.set_y(i, y);
    dtm.set_z(i, z);
  }

  return dtm;
}

}

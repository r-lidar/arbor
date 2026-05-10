/**
 * @file smooth3d.cpp
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
#include "nanoflann.h"

namespace arbor::utils {

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 3>;

PointCloud smooth3d(const PointCloud& cloud, double radius, int ncores)
{
  const size_t npoints = cloud.size();

  KDTree tree(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  tree.buildIndex();

  // Allocate output cloud (xyz only)
  PointCloud out(npoints, false);

  #pragma omp parallel for num_threads(ncores)
  for (int i = 0; i < static_cast<int>(npoints); i++)
  {
    double query[3];
    cloud.get_point(i, query);

    std::vector<nanoflann::ResultItem<uint32_t, double>> matches;
    nanoflann::SearchParameters params;
    params.sorted = false;

    tree.radiusSearch(query, radius, matches, params);

    out.set_x(i, cloud.get_x(i));
    out.set_y(i, cloud.get_y(i));
    out.set_z(i, cloud.get_z(i));

    double xtot = 0, ytot = 0, ztot = 0;
    for (const auto& m : matches)
    {
      const size_t idx = m.first;
      xtot += cloud.get_x(idx);
      ytot += cloud.get_y(idx);
      ztot += cloud.get_z(idx);
    }

    const double inv = 1.0 / static_cast<double>(matches.size());
    out.set_x(i, xtot * inv);
    out.set_y(i, ytot * inv);
    out.set_z(i, ztot * inv);
  }

  return out;
}

}

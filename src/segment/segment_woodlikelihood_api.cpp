/**
 * @file segment_woodlikelihood_api.cpp
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

#include <stdexcept>
#include <vector>

#include "arbor.h"
#include "PointCloud.h"

namespace arbor::segment {

void wood_likelihood(PointCloud& cloud, const settings::ArborParameters& params)
{
  if (!cloud.has_pwood()) throw std::runtime_error("wood_likelihood: point cloud is missing required 'pwood' attribute.");

  // The points below the cut height are classified foliage by the semantic segmentation whatever their wood
  // likelihood: they are hidden and get 0. The points within a small margin below the cut height are kept as
  // neighbours, so that the neighbourhoods of the points just above the cut are not truncated. Without
  // height above ground, every point is kept.
  constexpr double margin = 0.05;
  const double cut = params.global.cut_above_ground;
  const bool has_hag = cloud.has_hag();

  for (size_t i = 0 ; i < cloud.size() ; i++) cloud.set_pwood(i, 0);

  ScopedPartition<PointCloud> partition(cloud, [&](size_t i) { return !has_hag || cloud.get_hag(i) > cut - margin; });

  if (cloud.size() < static_cast<size_t>(params.woodlikelihood.k))
    throw std::runtime_error("wood_likelihood: fewer points above 'cut_above_ground' than 'k'");

  // pwood is indexed in the order of the partition: it is written in the point cloud within the scope to
  // follow the points back to their original index
  std::vector<float> pwood = utils::anisotropy(cloud, params.woodlikelihood.k);
  for (size_t i = 0 ; i < cloud.size() ; i++)
  {
    if (!has_hag || cloud.get_hag(i) > cut) cloud.set_pwood(i, pwood[i]);
  }
}

}

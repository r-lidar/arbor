/**
 * @file QSMbuilder_cluster.cpp
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

#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "arbor.h"
#include "QSMbuilder.h"
#include "dbscan/dbscan.hpp"

namespace arbor::qsm {

std::vector<int> QSMbuilder::cluster(const PointCloud& data, const std::vector<int>& iter, float eps)
{
  int minPts = 1;

  size_t n = data.size();
  std::vector<int> cl(n, 0);

  // Group point indices by their iter value
  std::map<int, std::vector<size_t>> groups;
  for (size_t i = 0; i < n; ++i)
  {
    groups[iter[i]].push_back(i);
  }

  // Unique pairs of <iter, dbscan_id>
  std::map<std::pair<int, int>, int> global_id_map;
  int next_global_id = 1;

  for (auto const& [group_id, indices] : groups)
  {
    // Prepare data for DBSCAN for this specific group
    std::vector<point3> points_in_group;
    points_in_group.reserve(indices.size());
    for (size_t idx : indices)
    {
      points_in_group.push_back({data.get_x(idx), data.get_y(idx), data.get_z(idx)});
    }

    // DBSCAN
    auto clusters = dbscan(points_in_group, eps, minPts);

    // Map relative indices back to local cluster labels
    // Local cluster 0 is noise
    std::vector<int> local_labels(indices.size(), 0);
    for (int c_idx = 0; c_idx < static_cast<int>(clusters.size()); ++c_idx)
    {
      int cluster_label = c_idx + 1; // Cluster IDs start at 1
      for (size_t relative_idx : clusters[c_idx])
      {
        local_labels[relative_idx] = cluster_label;
      }
    }

    // Assign global IDs based on (group_id, local_label)
    for (size_t i = 0; i < indices.size(); ++i)
    {
      std::pair<int, int> key = {group_id, local_labels[i]};

      if (global_id_map.find(key) == global_id_map.end())
      {
        global_id_map[key] = next_global_id++;
      }

      cl[indices[i]] = global_id_map[key];
    }
  }

  return cl;
}

}

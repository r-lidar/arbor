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

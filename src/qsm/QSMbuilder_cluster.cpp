#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "arbor.h"
#include "QSMbuilder.h"
#include "dbscan.hpp"

namespace arbor::qsm {

std::vector<std::pair<int, double>> QSMbuilder::clusters(const PointCloud& data, const std::vector<std::pair<int, double>>& iter_dist,  double cl_dist, const Logger& logger)
{
  logger("Computing clusters");

  // Result vectors for cluster ID and radius for each point
  std::vector<std::pair<int, double>> result(data.size(), {0, 0.0});

  // Get unique iterations and sort them
  std::set<int> unique_iters;
  for (const auto& [iter, dist] : iter_dist) {
    unique_iters.insert(iter);
  }

  bool first = true;
  double cl_d = cl_dist;

  for (int i : unique_iters)
  {
    // Find indices where iter == i
    std::vector<size_t> in_iter;
    for (size_t idx = 0; idx < iter_dist.size(); ++idx)
    {
      if (iter_dist[idx].first == i) {
        in_iter.push_back(idx);
      }
    }

    if (in_iter.empty()) continue;

    if (first)
    {
      // At the first iteration, only one cluster (the tree base)
      for (size_t idx : in_iter)
      {
        result[idx].first = 1;
      }

      // Compute the center of the first layer
      double mean_x = 0.0, mean_y = 0.0, mean_z = 0.0;
      for (size_t idx : in_iter)
      {
        mean_x += data.get_x(idx);
        mean_y += data.get_y(idx);
        mean_z += data.get_z(idx);
      }
      mean_x /= in_iter.size();
      mean_y /= in_iter.size();
      mean_z /= in_iter.size();

      // Compute radius for each point
      double radius_sum = 0.0;
      for (size_t idx : in_iter)
      {
        double dx = data.get_x(idx) - mean_x;
        double dy = data.get_y(idx) - mean_y;
        double dz = data.get_z(idx) - mean_z;
        double radius = std::sqrt(dx*dx + dy*dy + dz*dz);
        result[idx].second = radius;
        radius_sum += radius;
      }

      // The average radius is the first clustering distance
      cl_d = radius_sum / in_iter.size();
      first = false;
    }
    else
    {
      // Subsequent iterations
      if (in_iter.size() >= 2)
      {
        // Prepare data for dbscan
        std::vector<point3> points;
        points.reserve(in_iter.size());
        for (size_t idx : in_iter)
          points.push_back({data.get_x(idx), data.get_y(idx), data.get_z(idx)});

        // Perform clustering
        auto clusters = dbscan(points, cl_d, 1);

        // Assign cluster IDs
        for (size_t cluster_id = 0; cluster_id < clusters.size(); ++cluster_id)
        {
          for (size_t point_idx : clusters[cluster_id])
          {
            size_t original_idx = in_iter[point_idx];
            result[original_idx].first = static_cast<int>(cluster_id + 1);
          }
        }

        // Compute radius for each cluster
        std::map<int, std::vector<size_t>> cluster_map;
        for (size_t idx : in_iter)
        {
          int cluster_id = result[idx].first;
          cluster_map[cluster_id].push_back(idx);
        }

        double total_radius = 0.0;
        int radius_count = 0;

        for (const auto& [cluster_id, indices] : cluster_map)
        {
          // Compute center of this cluster
          double mean_x = 0.0, mean_y = 0.0, mean_z = 0.0;
          for (size_t idx : indices)
          {
            mean_x += data.get_x(idx);
            mean_y += data.get_y(idx);
            mean_z += data.get_z(idx);
          }
          mean_x /= indices.size();
          mean_y /= indices.size();
          mean_z /= indices.size();

          // Compute radius for each point in this cluster
          for (size_t idx : indices)
          {
            double dx = data.get_x(idx) - mean_x;
            double dy = data.get_y(idx) - mean_y;
            double dz = data.get_z(idx) - mean_z;
            double radius = std::sqrt(dx*dx + dy*dy + dz*dz);
            result[idx].second = radius;
            total_radius += radius;
            radius_count++;
          }
        }

        // Update clustering distance (average radius / 2)
        if (radius_count > 0)
        {
          //cl_d = (total_radius / radius_count) / 2.0;
        }
      }
      else
      {
        // Only one point, only one cluster
        result[in_iter[0]].first = 1;
        result[in_iter[0]].second = 0.0;
      }
    }
  }

  return result;
}

}

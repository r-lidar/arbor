#include "QSMbuilder.h"

#include <vector>
#include <cmath>
#include <array>
#include <unordered_map>
#include <limits>

std::vector<std::pair<int, double>> QSMbuilder::layers(const PointCloud& points, double D, const Logger& logger)
{
  struct GridKey {
    long x, y, z;
    bool operator==(const GridKey &other) const { return x == other.x && y == other.y && z == other.z; }
  };

  struct KeyHash {
    std::size_t operator()(const GridKey& k) const{
      return ((std::hash<long>()(k.x) ^ (std::hash<long>()(k.y) << 1)) >> 1) ^ (std::hash<long>()(k.z) << 1);
    }
  };

  logger("Computing layers");

  size_t n = points.size();
  double D2 = D * D;
  double invD = 1.0 / D; // Optimization: multiply is faster than divide

  // --- Spatial Hashing ---
  std::unordered_map<GridKey, std::vector<size_t>, KeyHash> grid;
  grid.reserve(n/100);

  double minZ = std::numeric_limits<double>::max();

  for (size_t i = 0 ; i < n ; i++)
  {
    double x = points.get_x(i);
    double y = points.get_y(i);
    double z = points.get_z(i);
    if (z < minZ) minZ = z;
    long gx = static_cast<long>(std::floor(x * invD));
    long gy = static_cast<long>(std::floor(y * invD));
    long gz = static_cast<long>(std::floor(z * invD));
    grid[{gx, gy, gz}].push_back(i);
  }

  // --- Output Initialization ---

  std::vector<std::pair<int, double>> result(n, {-1, std::numeric_limits<double>::quiet_NaN()});

  std::vector<size_t> current_layer;
  std::vector<size_t> next_layer;

  // Status tracker: -1 = unvisited, 0 = pending (in queue), 1 = processed
  std::vector<int8_t> status(n, -1);

  // --- Layer 1 Initialization ---
  // Select all points at the lowest Z level
  for (size_t i = 0 ; i < n ; i++)
  {
    double z = points.get_z(i);
    if (z <= minZ + 0.1)
    {
      result[i].first = 1;
      result[i].second = 0.0;
      status[i] = 1;
      current_layer.push_back(i);
    }
  }

  int current_iter_num = 2;
  size_t visited_count = current_layer.size();
  size_t search_start_idx = 0; // Optimization for restarting disconnected components

  // --- Propagation Loop ---
  while (visited_count < n)
  {
    // Handle Disconnected Components:
    // If the flood fill runs dry but points remain, pick the next unvisited point.
    if (current_layer.empty())
    {
      bool found_restart = false;
      for (size_t i = search_start_idx; i < n; ++i)
      {
        if (status[i] == -1)
        {
          result[i].first = current_iter_num;
          result[i].second = std::numeric_limits<double>::quiet_NaN(); // Start of new cluster
          status[i] = 1;
          current_layer.push_back(i);
          visited_count++;
          search_start_idx = i + 1; // Don't scan these again
          found_restart = true;
          break;
        }
      }
      if (!found_restart) break; // Should not happen given loop condition
    }

    next_layer.clear();

    // Process current wave
    for (size_t idx : current_layer)
    {
      double x = points.get_x(idx);
      double y = points.get_y(idx);
      double z = points.get_z(idx);

      long gx = static_cast<long>(std::floor(x * invD));
      long gy = static_cast<long>(std::floor(y * invD));
      long gz = static_cast<long>(std::floor(z * invD));

      // Check 3x3x3 Neighborhood
      for (long dx = -1; dx <= 1; ++dx)
      {
        for (long dy = -1; dy <= 1; ++dy)
        {
          for (long dz = -1; dz <= 1; ++dz)
          {
            GridKey key = {gx + dx, gy + dy, gz + dz};
            auto it = grid.find(key);
            if (it == grid.end()) continue;

            for (size_t neighbor_idx : it->second)
            {
              if (status[neighbor_idx] == 1) continue;

              double xn = points.get_x(neighbor_idx);
              double yn = points.get_y(neighbor_idx);
              double zn = points.get_z(neighbor_idx);
              double d2 = (x - xn)*(x - xn) + (y - yn)*(y - yn) + (z - zn)*(z - zn);

              if (d2 <= D2)
              {
                double d = std::sqrt(d2);

                if (status[neighbor_idx] == -1)
                {
                  // Found new point
                  status[neighbor_idx] = 0;
                  result[neighbor_idx].first = current_iter_num;
                  result[neighbor_idx].second = d;
                  next_layer.push_back(neighbor_idx);
                  visited_count++;
                }
                else if (status[neighbor_idx] == 0)
                {
                  // Point already in queue, check if this parent is closer
                  if (d < result[neighbor_idx].second)
                  {
                    result[neighbor_idx].second = d;
                  }
                }
              }
            }
          }
        }
      }
    }

    // Mark next layer as processed
    for (size_t idx : next_layer)
    {
      status[idx] = 1;
    }

    current_layer = std::move(next_layer);

    if (!current_layer.empty())
    {
      current_iter_num++;
    }
  }

  return result;
}

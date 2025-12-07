#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <array>
#include <unordered_map>
#include <limits>
#include <functional>

struct Point3D
{
  double x, y, z;
  size_t id;
};

struct GridKey
{
  long x, y, z;
  bool operator==(const GridKey &other) const
  {
    return x == other.x && y == other.y && z == other.z;
  }
};

struct KeyHash
{
  std::size_t operator()(const GridKey& k) const
  {
    return ((std::hash<long>()(k.x) ^ (std::hash<long>()(k.y) << 1)) >> 1) ^ (std::hash<long>()(k.z) << 1);
  }
};

struct LayerResult
{
  std::vector<int> iter;
  std::vector<double> dist;
};

LayerResult compute_layers(const std::vector<Point3D>& points, double D, std::function<void()> check_interrupt = nullptr)
{
  size_t n = points.size();
  double D2 = D * D;
  double invD = 1.0 / D; // Optimization: multiply is faster than divide

  // --- Spatial Hashing ---
  std::unordered_map<GridKey, std::vector<size_t>, KeyHash> grid;
  grid.reserve(n/100);

  double minZ = std::numeric_limits<double>::max();

  for (const auto& p : points)
  {
    if (p.z < minZ) minZ = p.z;
    long gx = static_cast<long>(std::floor(p.x * invD));
    long gy = static_cast<long>(std::floor(p.y * invD));
    long gz = static_cast<long>(std::floor(p.z * invD));
    grid[{gx, gy, gz}].push_back(p.id);
  }

  // --- Output Initialization ---

  LayerResult result;
  result.iter.assign(n, -1);
  result.dist.assign(n, std::numeric_limits<double>::quiet_NaN());

  std::vector<size_t> current_layer;
  std::vector<size_t> next_layer;

  // Status tracker: -1 = unvisited, 0 = pending (in queue), 1 = processed
  std::vector<int8_t> status(n, -1);

  // --- Layer 1 Initialization ---
  // Select all points at the lowest Z level
  for (const auto& p : points)
  {
    if (p.z <= minZ + 0.1)
    {
      result.iter[p.id] = 1;
      status[p.id] = 1;
      result.dist[p.id] = 0.0;
      current_layer.push_back(p.id);
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
          result.iter[i] = current_iter_num;
          status[i] = 1;
          result.dist[i] = std::numeric_limits<double>::quiet_NaN(); // Start of new cluster
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
      const Point3D& p = points[idx];

      long gx = static_cast<long>(std::floor(p.x * invD));
      long gy = static_cast<long>(std::floor(p.y * invD));
      long gz = static_cast<long>(std::floor(p.z * invD));

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
              if (status[neighbor_idx] == 1) continue; // Already done

              const Point3D& n_pt = points[neighbor_idx];
              double d2 = (p.x - n_pt.x)*(p.x - n_pt.x) +
                (p.y - n_pt.y)*(p.y - n_pt.y) +
                (p.z - n_pt.z)*(p.z - n_pt.z);

              if (d2 <= D2)
              {
                double d = std::sqrt(d2);

                if (status[neighbor_idx] == -1)
                {
                  // Found new point
                  status[neighbor_idx] = 0;
                  result.iter[neighbor_idx] = current_iter_num;
                  result.dist[neighbor_idx] = d;
                  next_layer.push_back(neighbor_idx);
                  visited_count++;
                }
                else if (status[neighbor_idx] == 0)
                {
                  // Point already in queue, check if this parent is closer
                  if (d < result.dist[neighbor_idx])
                  {
                    result.dist[neighbor_idx] = d;
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

    if (check_interrupt) check_interrupt();
  }

  return result;
}

// -------------------------------------------------------------------------
// RCPP INTERFACE SECTION
// -------------------------------------------------------------------------

Rcpp::DataFrame qsm_layers_cpp(Rcpp::DataFrame df, double D)
{
  // 1. Unpack R DataFrame
  Rcpp::NumericVector X = df["X"];
  Rcpp::NumericVector Y = df["Y"];
  Rcpp::NumericVector Z = df["Z"];

  int n = X.size();

  // 2. Convert to Pure C++ Vector
  std::vector<Point3D> points(n);
  for(int i = 0; i < n; ++i)
  {
    points[i] = { X[i], Y[i], Z[i], (size_t)i };
  }

  // 3. Call Core Logic
  // We pass Rcpp::checkUserInterrupt as a lambda to keep core logic pure
  auto interrupt_callback = []() { Rcpp::checkUserInterrupt(); };
  LayerResult result = compute_layers(points, D, interrupt_callback);

  // 4. Repack into R DataFrame
  // We reuse the input columns to save memory/time, or create new ones if needed.
  // Here we create a new DataFrame as requested.
  return Rcpp::DataFrame::create(
    Rcpp::_["X"] = X,
    Rcpp::_["Y"] = Y,
    Rcpp::_["Z"] = Z,
    Rcpp::_["iter"] = result.iter,
    Rcpp::_["dist"] = result.dist
  );
}



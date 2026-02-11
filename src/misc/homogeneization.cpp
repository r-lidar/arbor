#include <unordered_map>
#include <vector>
#include <random>

#include "Adaptor.h"

std::vector<bool> homogeneization(const PointCloud& pc, double res, bool hybrid = true)
{

  size_t n = pc.n_points;

  // COMPUTE VOXEL ID FOR EACH POINT
  // ===============================

  // --- Compute bounding box ---
  double xmin = pc.get_x(0), xmax = pc.get_x(0);
  double ymin = pc.get_y(0), ymax = pc.get_y(0);
  double zmin = pc.get_z(0), zmax = pc.get_z(0);

  for (size_t i = 1; i < n; ++i)
  {
    double x = pc.get_x(i), y = pc.get_y(i), z = pc.get_z(i);
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
  }

  // --- Compute grid dimensions ---
  int64_t nx = static_cast<int64_t>(std::floor((xmax - xmin) / res)) + 1;
  int64_t ny = static_cast<int64_t>(std::floor((ymax - ymin) / res)) + 1;
  int64_t nz = static_cast<int64_t>(std::floor((zmax - zmin) / res)) + 1;

  // --- Compute voxel IDs ---
  std::vector<int64_t> id(n);
  for (size_t i = 0; i < n; ++i)
  {
    int64_t ix = static_cast<int64_t>(std::floor((pc.get_x(i) - xmin) / res));
    int64_t iy = static_cast<int64_t>(std::floor((pc.get_y(i) - ymin) / res));
    int64_t iz = static_cast<int64_t>(std::floor((pc.get_z(i) - zmin) / res));
    id[i] = ((iz * ny) + iy) * nx + ix;
  }

  // --- Deduplicate IDs ---
  std::vector<int64_t> sorted = id;
  std::sort(sorted.begin(), sorted.end());
  sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());

  // --- Remap to dense integers using binary search ---
  std::vector<int64_t> remapped(n);
  for (size_t i = 0; i < n; ++i)
  {
    int64_t idx = std::lower_bound(sorted.begin(), sorted.end(), id[i]) - sorted.begin();
    remapped[i] = idx;
  }

  std::swap(id, remapped);

  // SELECT POINTS TO RETAIN
  // (Barycentric decimation)
  // ===============================

  // Step 1: group indices by voxel ID
  std::unordered_map<int, std::vector<int>> groups;
  groups.reserve(n / 10);

  for (int i = 0; i < n; ++i)
  {
    int key = (int)id[i];
    groups[key].push_back(i);
  }

  // Step 2: result vector (initialized to FALSE)
  std::vector<bool> keep(n, false);

  // Step 3: for each voxel group, find closest to mean
  for (auto &kv : groups)
  {
    const std::vector<int> &idx = kv.second;
    const int m = idx.size();

    if (m == 1)
    {
      keep[idx[0]] = true;
      continue;
    }

    // compute mean
    double mx = 0, my = 0, mz = 0;
    for (int j : idx)
    {
      mx += pc.get_x(j);
      my += pc.get_y(j);
      mz += pc.get_z(j);
    }
    mx /= m;
    my /= m;
    mz /= m;

    // find closest
    double bestDist = R_PosInf;
    int bestIdx = idx[0];
    for (int j : idx)
    {
      double dx = pc.get_x(j) - mx;
      double dy = pc.get_y(j) - my;
      double dz = pc.get_z(j) - mz;
      double d = dx * dx + dy * dy + dz * dz;
      if (d < bestDist) {
        bestDist = d;
        bestIdx = j;
      }
    }

    keep[bestIdx] = true;
  }

  if (!hybrid) return keep;

  // ===============================
  // (HYBRID DECIMATION
  // Random reinjection of non sampled point
  // ===============================

  // Count TRUEs
  n = std::count(keep.begin(), keep.end(), true);

  // Collect indices where keep == false
  std::vector<size_t> rm;
  rm.reserve(keep.size());
  for (size_t i = 0; i < keep.size(); ++i)
  {
    if (!keep[i]) {
      rm.push_back(i);
    }
  }

  // Number of indices to sample
  size_t k = std::min(static_cast<size_t>(n * 0.1), rm.size());
  if (k == 0) return keep;

  // Random generator
  static std::mt19937 rng(std::random_device{}());

  // Shuffle and take first k
  std::shuffle(rm.begin(), rm.end(), rng);

  for (size_t i = 0; i < k; ++i)
  {
    keep[rm[i]] = true;
  }

  return keep;
}

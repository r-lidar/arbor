#include <unordered_map>
#include <vector>
#include <random>
#include <algorithm>

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

  // OPTIMIZATION 1: Use counting sort instead of std::sort for better performance
  // when the number of unique voxels is much smaller than n
  // BUT: Only if we expect many duplicates. Otherwise fallback to original.
  // Let's optimize the remap step instead.

  // OPTIMIZATION 2: Combined sort and remap using indices
  // Create index array to avoid multiple passes
  std::vector<int> indices(n);
  for (size_t i = 0; i < n; ++i) indices[i] = i;

  // Sort indices by voxel ID (not the IDs themselves)
  std::sort(indices.begin(), indices.end(), [&id](int a, int b) {
    return id[a] < id[b];
  });

  // Now points with same voxel ID are adjacent in the sorted order
  // We can build groups AND do remapping in one pass

  std::vector<bool> keep(n, false);

  size_t i = 0;
  while (i < n)
  {
    int64_t current_voxel = id[indices[i]];
    size_t group_start = i;

    // Find end of this voxel group
    while (i < n && id[indices[i]] == current_voxel) {
      ++i;
    }

    size_t group_size = i - group_start;

    if (group_size == 1)
    {
      keep[indices[group_start]] = true;
      continue;
    }

    // Compute mean for this group
    double mx = 0, my = 0, mz = 0;
    for (size_t j = group_start; j < i; ++j)
    {
      int idx = indices[j];
      mx += pc.get_x(idx);
      my += pc.get_y(idx);
      mz += pc.get_z(idx);
    }
    mx /= group_size;
    my /= group_size;
    mz /= group_size;

    // Find closest point to mean
    double bestDist = R_PosInf;
    int bestIdx = indices[group_start];
    for (size_t j = group_start; j < i; ++j)
    {
      int idx = indices[j];
      double dx = pc.get_x(idx) - mx;
      double dy = pc.get_y(idx) - my;
      double dz = pc.get_z(idx) - mz;
      double d = dx * dx + dy * dy + dz * dz;
      if (d < bestDist) {
        bestDist = d;
        bestIdx = idx;
      }
    }

    keep[bestIdx] = true;
  }

  if (!hybrid) return keep;

  // ===============================
  // HYBRID DECIMATION
  // ===============================

  size_t kept = std::count(keep.begin(), keep.end(), true);

  std::vector<size_t> rm;
  rm.reserve(n - kept);
  for (size_t i = 0; i < n; ++i)
  {
    if (!keep[i]) {
      rm.push_back(i);
    }
  }

  size_t k = std::min(static_cast<size_t>(kept * 0.1), rm.size());
  if (k == 0) return keep;

  static std::mt19937 rng(std::random_device{}());
  std::shuffle(rm.begin(), rm.end(), rng);

  for (size_t i = 0; i < k; ++i)
  {
    keep[rm[i]] = true;
  }

  return keep;
}

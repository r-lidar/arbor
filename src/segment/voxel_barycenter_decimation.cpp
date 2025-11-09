#include <Rcpp.h>
#include <unordered_map>
#include <vector>
using namespace Rcpp;


LogicalVector C_voxel_barycenter_decimate(NumericVector X, NumericVector Y,  NumericVector Z, double res)
{

  size_t n = X.size();
  if (Y.size() != n || Z.size() != n)
    stop("Input vectors must have same length");

  // COMPUTE VOXEL ID FOR EACH POINT
  // ===============================

  // --- Compute bounding box ---
  double xmin = X[0], xmax = X[0];
  double ymin = Y[0], ymax = Y[0];
  double zmin = Z[0], zmax = Z[0];

  for (size_t i = 1; i < n; ++i)
  {
    double x = X[i], y = Y[i], z = Z[i];
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
    int64_t ix = static_cast<int64_t>(std::floor((X[i] - xmin) / res));
    int64_t iy = static_cast<int64_t>(std::floor((Y[i] - ymin) / res));
    int64_t iz = static_cast<int64_t>(std::floor((Z[i] - zmin) / res));
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
  LogicalVector keep(n, false);

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
      mx += X[j];
      my += Y[j];
      mz += Z[j];
    }
    mx /= m;
    my /= m;
    mz /= m;

    // find closest
    double bestDist = R_PosInf;
    int bestIdx = idx[0];
    for (int j : idx)
    {
      double dx = X[j] - mx;
      double dy = Y[j] - my;
      double dz = Z[j] - mz;
      double d = dx * dx + dy * dy + dz * dz;
      if (d < bestDist) {
        bestDist = d;
        bestIdx = j;
      }
    }

    keep[bestIdx] = true;
  }

  return keep;
}

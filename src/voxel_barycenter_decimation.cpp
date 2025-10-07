#include <Rcpp.h>
#include <unordered_map>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector C_voxel_barycenter_decimate(
    NumericVector X,
    NumericVector Y,
    NumericVector Z,
    NumericVector id)
{
  const int n = X.size();
  if (Y.size() != n || Z.size() != n || id.size() != n)
    stop("Input vectors must have same length");

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

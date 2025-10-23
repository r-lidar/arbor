#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <algorithm>
#include <sstream>

using namespace Rcpp;

struct ClusterCenter
{
  double x, y, z;
  int iter;
  int id;
  bool done = false;
};

struct pair_hash
{
  std::size_t operator()(const std::pair<int, int>& p) const
  {
    return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
  }
};

DataFrame cpp_build_skeleton_old(DataFrame data, double max_d)
{
  NumericVector X = data["X"];
  NumericVector Y = data["Y"];
  NumericVector Z = data["Z"];
  IntegerVector iter = data["iter"];
  IntegerVector cluster = data["cluster"];

  // Step 1: compute mean for each (iter, cluster)
  typedef std::pair<int, int> ClusterKey;
  struct ClusterSum
  {
    double x = 0, y = 0, z = 0;
    int count = 0;
  };

  std::unordered_map<ClusterKey, ClusterSum, pair_hash> cluster_sums;
  for (int i = 0; i < X.size(); ++i)
  {
    ClusterKey key = std::make_pair(iter[i], cluster[i]);
    cluster_sums[key].x += X[i];
    cluster_sums[key].y += Y[i];
    cluster_sums[key].z += Z[i];
    cluster_sums[key].count++;
  }

  // Build cluster centers
  std::vector<ClusterCenter> centers;
  std::unordered_map<int, ClusterCenter*> centerByID;
  int id = 1;
  for (auto& entry : cluster_sums)
  {
    auto& key = entry.first;
    auto& sum = entry.second;
    ClusterCenter center;
    center.x = sum.x / sum.count;
    center.y = sum.y / sum.count;
    center.z = sum.z / sum.count;
    center.iter = key.first;
    center.id = id++;
    centers.push_back(center);
  }

  // Copy centers for the search set
  std::vector<ClusterCenter*> searchSpace;
  for (auto& c : centers)
  {
    centerByID[c.id] = &c;
    searchSpace.push_back(&c);
  }

  // Find initial root: min Z
  ClusterCenter* root = *std::min_element(searchSpace.begin(), searchSpace.end(), [](ClusterCenter* a, ClusterCenter* b) { return a->z < b->z; });
  root->done = true;
  searchSpace.erase(std::remove(searchSpace.begin(), searchSpace.end(), root), searchSpace.end());

  std::vector<double> startX, startY, startZ, endX, endY, endZ;
  const double max_d2 = max_d * max_d;

  while (!searchSpace.empty())
  {
    ClusterCenter* start = root;
    ClusterCenter* newRoot = nullptr;
    double bestD2 = std::numeric_limits<double>::max();

    for (auto* c : searchSpace)
    {
      if (c->iter <= root->iter) continue;
      double dx = c->x - root->x;
      double dy = c->y - root->y;
      double dz = c->z - root->z;
      double d2 = dx*dx + dy*dy + dz*dz;
      if (d2 < bestD2 && d2 <= max_d2)
      {
        bestD2 = d2;
        newRoot = c;
      }
    }

    if (newRoot)
    {
      newRoot->done = true;
      searchSpace.erase(std::remove(searchSpace.begin(), searchSpace.end(), newRoot), searchSpace.end());
      startX.push_back(start->x);
      startY.push_back(start->y);
      startZ.push_back(start->z);
      endX.push_back(newRoot->x);
      endY.push_back(newRoot->y);
      endZ.push_back(newRoot->z);
      root = newRoot;
    }
    else
    {
      // fallback: take lowest iter
      auto minIt = std::min_element(searchSpace.begin(), searchSpace.end(), [](ClusterCenter* a, ClusterCenter* b) { return a->iter < b->iter; });
      root = *minIt;

      // find nearest "done"
      ClusterCenter* nearestDone = nullptr;
      double bestDist = std::numeric_limits<double>::max();
      for (auto& c : centers)
      {
        if (!c.done) continue;
        double dx = c.x - root->x;
        double dy = c.y - root->y;
        double dz = c.z - root->z;
        double d = dx*dx + dy*dy + dz*dz;
        if (d < bestDist)
        {
          bestDist = d;
          nearestDone = &c;
        }
      }

      if (!nearestDone) break;

      root->done = true;
      searchSpace.erase(minIt);
      startX.push_back(nearestDone->x);
      startY.push_back(nearestDone->y);
      startZ.push_back(nearestDone->z);
      endX.push_back(root->x);
      endY.push_back(root->y);
      endZ.push_back(root->z);
    }
  }

  return DataFrame::create(
    _["startX"] = startX,
    _["startY"] = startY,
    _["startZ"] = startZ,
    _["endX"] = endX,
    _["endY"] = endY,
    _["endZ"] = endZ
  );
}

#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <algorithm>
#include <sstream>

#include "ransac.h"

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

DataFrame cpp_build_skeleton(DataFrame data, double max_d)
{
  // Extract numeric columns from the input R DataFrame.
  // Each vector corresponds to a point in 3D space, associated with an iteration and a cluster ID.
  NumericVector X = data["X"];
  NumericVector Y = data["Y"];
  NumericVector Z = data["Z"];
  IntegerVector iter = data["iter"];
  IntegerVector cluster = data["cluster"];

  // Step 1: compute mean for each (iter, cluster)
  // The code groups all points sharing the same combination of (iter, cluster)
  // and computes their mean X, Y, Z coordinates — i.e., the cluster "centroid" for that iteration.

  typedef std::pair<int, int> ClusterKey;
  std::unordered_map<ClusterKey, std::vector<int>, pair_hash> cluster_indices;

  // group indices by (iter, cluster)
  for (int i = 0; i < X.size(); ++i)
  {
    cluster_indices[std::make_pair(iter[i], cluster[i])].push_back(i);
  }

  // build centers using RANSAC and our new class
  std::vector<ClusterCenter> centers;
  int id = 1;

  for (auto& entry : cluster_indices)
  {
    auto key = entry.first;
    auto& indices = entry.second;

    ClusterCenter c;
    c.iter = key.first;
    c.id = id++;

    if (indices.size() >= 100)
    {
      RansacCircle rc(100, 0.02); // iterations, inlier threshold, early exit ratio

      for (int idx : indices) rc.add_point(X[idx], Y[idx], Z[idx]);

      rc.find_circle();

      if (rc.is_valid(0.5, 0.3, 120.0))
      {
        auto center = rc.get_center();
        c.x = center[0];
        c.y = center[1];
        c.z = center[2];
      }
      else
      {
        // fallback to average if RANSAC fails validation
        double sumx = 0, sumy = 0, sumz = 0;
        for (int idx : indices)
        {
          sumx += X[idx];
          sumy += Y[idx];
          sumz += Z[idx];
        }
        c.x = sumx / indices.size();
        c.y = sumy / indices.size();
        c.z = sumz / indices.size();
      }
    }
    else
    {
      // fallback to average for small clusters
      double sumx = 0, sumy = 0, sumz = 0;
      for (int idx : indices)
      {
        sumx += X[idx];
        sumy += Y[idx];
        sumz += Z[idx];
      }
      c.x = sumx / indices.size();
      c.y = sumy / indices.size();
      c.z = sumz / indices.size();
    }

    centers.push_back(c);
  }

  // Prepare for neighbor searching.
  // We create a vector of pointers to the actual centers (so we can modify them in place).
  std::vector<ClusterCenter*> searchSpace;
  for (auto& c : centers)
    searchSpace.push_back(&c);

  // Find initial root: cluster center with minimum Z value (lowest in space).
  // This acts as the starting point for building the skeleton structure.
  ClusterCenter* root = *std::min_element(
    searchSpace.begin(),
    searchSpace.end(),
    [](ClusterCenter* a, ClusterCenter* b) { return a->z < b->z; }
  );

  root->done = true;  // Mark root as processed
  // Remove the root from the search pool (it’s now part of the skeleton).
  searchSpace.erase(std::remove(searchSpace.begin(), searchSpace.end(), root), searchSpace.end());

  // Prepare vectors to store edges of the skeleton.
  // Each edge connects one cluster center (start) to another (end).
  std::vector<double> startX, startY, startZ, endX, endY, endZ;

  // Precompute max distance squared for faster comparison.
  const double max_d2 = max_d * max_d;

  // Main loop: keep connecting cluster centers until all are processed.
  while (!searchSpace.empty())
  {
    ClusterCenter* start = root;          // Current "root" center (start point)
    ClusterCenter* newRoot = nullptr;     // Will hold the nearest valid neighbor
    double bestD2 = std::numeric_limits<double>::max(); // Track the smallest distance squared found

    // Search among remaining centers for the next connection
    for (auto* c : searchSpace)
    {
      // Only connect centers that belong to a later iteration (enforces directionality)
      if (c->iter <= root->iter) continue;

      // Compute squared Euclidean distance between centers
      double dx = c->x - root->x;
      double dy = c->y - root->y;
      double dz = c->z - root->z;
      double d2 = dx*dx + dy*dy + dz*dz;

      // Keep the closest candidate under the distance threshold
      if (d2 < bestD2 && d2 <= max_d2)
      {
        bestD2 = d2;
        newRoot = c;  // Best candidate so far
      }
    }

    if (newRoot)
    {
      // Found a valid neighbor close enough — connect root → newRoot
      newRoot->done = true;  // Mark as processed
      // Remove it from the search pool
      searchSpace.erase(std::remove(searchSpace.begin(), searchSpace.end(), newRoot), searchSpace.end());

      // Record the edge from root to newRoot
      startX.push_back(start->x);
      startY.push_back(start->y);
      startZ.push_back(start->z);
      endX.push_back(newRoot->x);
      endY.push_back(newRoot->y);
      endZ.push_back(newRoot->z);

      // Move to the next root (traverse forward)
      root = newRoot;
    }
    else
    {
      // No nearby candidate found within distance threshold.
      // Fallback strategy: pick the point with the lowest iteration value among unprocessed ones.
      auto minIt = std::min_element(
        searchSpace.begin(),
        searchSpace.end(),
        [](ClusterCenter* a, ClusterCenter* b) { return a->iter < b->iter; }
      );
      root = *minIt;

      // Now find the nearest already-processed ("done") node in Euclidean space.
      ClusterCenter* nearestDone = nullptr;
      double bestDist = std::numeric_limits<double>::max();
      for (auto& c : centers)
      {
        if (!c.done) continue;  // Skip unprocessed centers

        // Distance to this processed center
        double dx = c.x - root->x;
        double dy = c.y - root->y;
        double dz = c.z - root->z;
        double d = dx*dx + dy*dy + dz*dz;

        if (d < bestDist)
        {
          bestDist = d;
          nearestDone = &c;  // Best processed neighbor
        }
      }

      // If no processed center exists, the skeleton cannot continue
      if (!nearestDone) break;

      // Connect new root to its nearest processed neighbor
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

  // Finally, build an R DataFrame of all the edges forming the skeleton.
  // Each row represents a connection between two cluster centers.
  return DataFrame::create(
    _["startX"] = startX,
    _["startY"] = startY,
    _["startZ"] = startZ,
    _["endX"] = endX,
    _["endY"] = endY,
    _["endZ"] = endZ
  );
}

#include <Rcpp.h>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <limits>
#include "nanoflann/nanoflann.h"

using namespace Rcpp;
using namespace nanoflann;

// Compact point cloud structure for nanoflann
struct SimpleCloud
{
  std::vector<std::array<double, 3>> pts;

  inline size_t kdtree_get_point_count() const { return pts.size(); }

  inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
    return pts[idx][dim];
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const { return false; }
};

DataFrame cpp_compute_layers(NumericMatrix coords, double D)
{
  const int n = coords.nrow();
  const double D2 = D * D; // precompute squared distance threshold

  IntegerVector ID(n), iter(n, -1);
  NumericVector dist(n, NA_REAL);

  for (int i = 0; i < n; ++i) ID[i] = i;

  // Identify the minimum Z value
  double minZ = coords(0, 2);
  for (int i = 1; i < n; ++i)
  {
    if (coords(i, 2) < minZ)
      minZ = coords(i, 2);
  }

  std::vector<int> layer;
  std::unordered_set<int> remaining;

  // Assign first layer
  for (int i = 0; i < n; ++i)
  {
    if (coords(i, 2) <= minZ + 0.1)
    {
      iter[i] = 1;
      layer.push_back(i);
    }
    else
    {
      remaining.insert(i);
    }
  }

  int current_iter = 2;

  while (!remaining.empty())
  {
    // Build point cloud for current layer
    SimpleCloud layer_cloud;
    layer_cloud.pts.reserve(layer.size());
    for (int idx : layer)
    {
      layer_cloud.pts.push_back({coords(idx, 0), coords(idx, 1), coords(idx, 2)});
    }

    typedef KDTreeSingleIndexAdaptor<
      L2_Simple_Adaptor<double, SimpleCloud>,
      SimpleCloud, 3
    > KDTree;

    // Construct KD-tree from current layer
    KDTree index(3, layer_cloud, KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    std::vector<int> next_layer, to_remove;

    // Find neighbors in current layer
    for (int idx : remaining)
    {
      double query_pt[3] = { coords(idx, 0), coords(idx, 1), coords(idx, 2) };

      size_t ret_index;
      double out_dist_sqr;
      nanoflann::KNNResultSet<double> resultSet(1);
      resultSet.init(&ret_index, &out_dist_sqr);
      index.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(10));

      dist[idx] = std::sqrt(out_dist_sqr);
      if (out_dist_sqr <= D2)
      {
        next_layer.push_back(idx);
        to_remove.push_back(idx);
      }
    }

    // Handle disconnected points by choosing closest
    if (next_layer.empty())
    {
      SimpleCloud ref_cloud;
      for (int i = 0; i < n; ++i)
      {
        if (iter[i] != -1)
          ref_cloud.pts.push_back({coords(i, 0), coords(i, 1), coords(i, 2)});
      }

      KDTree ref_index(3, ref_cloud, KDTreeSingleIndexAdaptorParams(10));
      ref_index.buildIndex();

      double min_dist = std::numeric_limits<double>::max();
      int closest_idx = -1;

      for (int idx : remaining)
      {
        double query_pt[3] = { coords(idx, 0), coords(idx, 1), coords(idx, 2) };
        size_t ret_index;
        double out_dist_sqr;
        nanoflann::KNNResultSet<double> resultSet(1);
        resultSet.init(&ret_index, &out_dist_sqr);
        ref_index.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(10));

        if (out_dist_sqr < min_dist)
        {
          min_dist = out_dist_sqr;
          closest_idx = idx;
        }
      }

      next_layer.push_back(closest_idx);
      to_remove.push_back(closest_idx);
      dist[closest_idx] = std::sqrt(min_dist);
    }

    // Update iteration index and remove processed points
    for (int idx : next_layer) iter[idx] = current_iter;
    for (int idx : to_remove) remaining.erase(idx);

    layer = std::move(next_layer);
    ++current_iter;
  }

  // Return result as DataFrame
  return DataFrame::create(
    _["X"] = coords(_, 0),
    _["Y"] = coords(_, 1),
    _["Z"] = coords(_, 2),
    _["ID"] = ID,
    _["iter"] = iter,
    _["dist"] = dist
  );
}

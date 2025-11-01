#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include "nanoflann/nanoflann.h"

using namespace Rcpp;
using namespace nanoflann;

// Compact point cloud for nanoflann
struct SimpleCloud {
  std::vector<std::array<double, 3>> pts;
  inline size_t kdtree_get_point_count() const { return pts.size(); }
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
    return pts[idx][dim];
  }
  template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
};

// [[Rcpp::export]]
DataFrame cpp_compute_layers(NumericMatrix coords, double D)
{
  const int n = coords.nrow();
  const double D2 = D * D;

  // Convert to efficient native format
  std::vector<std::array<double, 3>> points(n);
  for (int i = 0; i < n; ++i)
    points[i] = {coords(i,0), coords(i,1), coords(i,2)};

  IntegerVector ID(n), iter(n, -1);
  NumericVector dist(n, NA_REAL);
  for (int i = 0; i < n; ++i) ID[i] = i;

  // Find min Z
  double minZ = points[0][2];
  for (int i = 1; i < n; ++i)
    if (points[i][2] < minZ) minZ = points[i][2];

    std::vector<int> layer;
    std::vector<char> active(n, 1);  // active == 1 means still remaining

    // First layer (lowest Z)
    for (int i = 0; i < n; ++i) {
      if (points[i][2] <= minZ + 0.1) {
        iter[i] = 1;
        active[i] = 0;
        layer.push_back(i);
      }
    }

    int current_iter = 2;
    std::vector<int> next_layer;
    next_layer.reserve(n);

    SimpleCloud ref_cloud{points};
    typedef KDTreeSingleIndexAdaptor<
      L2_Simple_Adaptor<double, SimpleCloud>, SimpleCloud, 3> KDTree;

    // Single global KD-tree (reused for distance fallback)
    KDTree ref_index(3, ref_cloud, KDTreeSingleIndexAdaptorParams(10));
    ref_index.buildIndex();

    while (true) {
      // Build KD-tree for current layer
      SimpleCloud layer_cloud;
      layer_cloud.pts.reserve(layer.size());
      for (int idx : layer)
        layer_cloud.pts.push_back(points[idx]);

      KDTree index(3, layer_cloud, KDTreeSingleIndexAdaptorParams(10));
      index.buildIndex();

      next_layer.clear();
      bool any_active = false;

      for (int i = 0; i < n; ++i) {
        if (!active[i]) continue;

        double query_pt[3] = {points[i][0], points[i][1], points[i][2]};
        size_t ret_index;
        double out_dist_sqr;
        nanoflann::KNNResultSet<double> resultSet(1);
        resultSet.init(&ret_index, &out_dist_sqr);
        index.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters());

        if (out_dist_sqr <= D2) {
          iter[i] = current_iter;
          dist[i] = std::sqrt(out_dist_sqr);
          active[i] = 0;
          next_layer.push_back(i);
          any_active = true;
        }
      }

      // handle disconnected case
      if (next_layer.empty()) {
        int closest_idx = -1;
        double min_dist = std::numeric_limits<double>::max();

        for (int i = 0; i < n; ++i) {
          if (!active[i]) continue;
          double query_pt[3] = {points[i][0], points[i][1], points[i][2]};
          size_t ret_index;
          double out_dist_sqr;
          nanoflann::KNNResultSet<double> resultSet(1);
          resultSet.init(&ret_index, &out_dist_sqr);
          ref_index.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(10));
          if (out_dist_sqr < min_dist) {
            min_dist = out_dist_sqr;
            closest_idx = i;
          }
        }

        if (closest_idx == -1) break; // no remaining points
        iter[closest_idx] = current_iter;
        dist[closest_idx] = std::sqrt(min_dist);
        active[closest_idx] = 0;
        next_layer.push_back(closest_idx);
        any_active = true;
      }

      if (!any_active) break;
      layer.swap(next_layer);
      ++current_iter;
    }

    return DataFrame::create(
      _["X"] = coords(_, 0),
      _["Y"] = coords(_, 1),
      _["Z"] = coords(_, 2),
      _["ID"] = ID,
      _["iter"] = iter,
      _["dist"] = dist
    );
}

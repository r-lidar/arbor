#include <vector>
#include <atomic>
#include <cmath>

#include "myomp.h"
#include "nanoflann.h"
#include "Adaptor.h"
#include "progressbar.h"

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>,PointCloud, 3>;
using index_t = nanoflann::KNNResultSet<double>::IndexType;

std::vector<bool> sor(const PointCloud& pc, unsigned int k, double m, int ncpu)
{
  size_t n = pc.n_points;

  if (k <= 2) throw std::invalid_argument("k must be > 2");
  if (k >= n) throw std::invalid_argument("k must be < number of points");

  KDTree tree(3, pc, nanoflann::KDTreeSingleIndexAdaptorParams(40));
  tree.buildIndex();
  nanoflann::SearchParameters params;
  params.sorted = false;

  Progress pb(n, "SOR");
  std::atomic<bool> abort(false);
  std::vector<double> dmean(n);

  #pragma omp parallel num_threads(ncpu)
  {
    std::vector<index_t> idx(k);
    std::vector<double> dist(k);
    nanoflann::KNNResultSet<double> resultSet(k);
    double q[3];

    #pragma omp for schedule(static)
    for (size_t i = 0; i < n; ++i)
    {
      if (abort.load(std::memory_order_relaxed)) continue;
      if(pb.check_interrupt()) abort = true;
      pb.tick();

      q[0] = pc.get_x(i);
      q[1] = pc.get_y(i);
      q[2] = pc.get_z(i);

      resultSet.init(idx.data(), dist.data());
      tree.findNeighbors(resultSet, q, params);

      double dsum = 0.0;
      unsigned int cnt = 0;
      for (unsigned int j = 0; j < k; ++j)
      {
        if (idx[j] == i) continue;   // safest
        double d = std::sqrt(dist[j]);
        if (d > 0)
        {
          dsum += d;
          cnt++;
        }
      }

      dmean[i] = (cnt > 0) ? dsum / cnt : 0.0;
    }

    pb.tick();
  }

  if (abort.load()) Rcpp::stop("Computation aborted");

  double mean = 0.0;
  for (double v : dmean) mean += v;
  mean /= static_cast<double>(n);

  double var = 0.0;
  for (double v : dmean)
  {
    double dv = v - mean;
    var += dv * dv;
  }
  var /= static_cast<double>(n - 1);
  double sd = std::sqrt(var);

  double threshold = mean + m * sd;

  std::vector<bool> out(n);
  for (size_t i = 0; i < n; ++i)
    out[i] = (dmean[i] > threshold);

  return out;
}

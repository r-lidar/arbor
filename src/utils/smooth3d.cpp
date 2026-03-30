#include "arbor.h"
#include "nanoflann.h"

namespace arbor::utils {

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 3>;

PointCloud smooth3d(const PointCloud& cloud, double radius, int ncores)
{
  const size_t npoints = cloud.size();

  KDTree tree(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  tree.buildIndex();

  // Allocate output cloud (xyz only)
  PointCloud out(npoints, false);

  #pragma omp parallel for num_threads(ncores)
  for (int i = 0; i < static_cast<int>(npoints); i++)
  {
    double query[3];
    cloud.get_point(i, query);

    std::vector<nanoflann::ResultItem<uint32_t, double>> matches;
    nanoflann::SearchParameters params;
    params.sorted = false;

    tree.radiusSearch(query, radius, matches, params);

    out.set_x(i, cloud.get_x(i));
    out.set_y(i, cloud.get_y(i));
    out.set_z(i, cloud.get_z(i));

    double xtot = 0, ytot = 0, ztot = 0;
    for (const auto& m : matches)
    {
      const size_t idx = m.first;
      xtot += cloud.get_x(idx);
      ytot += cloud.get_y(idx);
      ztot += cloud.get_z(idx);
    }

    const double inv = 1.0 / static_cast<double>(matches.size());
    out.set_x(i, xtot * inv);
    out.set_y(i, ytot * inv);
    out.set_z(i, ztot * inv);
  }

  return out;
}

}

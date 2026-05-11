#include "QSMbuilder.h"
#include "PointCloud.h"
#include "nanoflann.h"

namespace arbor::qsm {

class SimpleAdaptor
{
public:
  struct Point { double x, y, z; int id; };
  std::vector<Point> points;
  inline size_t kdtree_get_point_count() const { return points.size(); }
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const
  {
    if (dim == 0) return points[idx].x;
    if (dim == 1) return points[idx].y;
    return points[idx].z;
  }
  template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
  inline size_t point_count() const { return points.size(); }
  inline size_t size() const { return points.size(); }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, SimpleAdaptor>, SimpleAdaptor, 3> CentroidKDTree;

std::unordered_map<int, std::vector<size_t>> QSMbuilder::group_points_by_edge(const QSM& graph, const PointCloud& tree)
{
  std::unordered_map<int, std::vector<size_t>> points_per_eid;
  if (graph.edge_count() == 0) return points_per_eid;

  // Prepare centroids for KD-Tree
  SimpleAdaptor centroids_cloud;
  centroids_cloud.points.reserve(graph.edge_count());

  // We need this to map the KD-Tree result index back to the actual Edge ID
  std::vector<int> index_to_eid;
  index_to_eid.reserve(graph.edge_count());

  for (const auto& [eid, einfo] : graph.edges())
  {
    const auto& src = graph.node(einfo.source);
    const auto& tgt = graph.node(einfo.target);

    double cx = (src.x + tgt.x) / 2.0;
    double cy = (src.y + tgt.y) / 2.0;
    double cz = (src.z + tgt.z) / 2.0;

    centroids_cloud.points.push_back({cx, cy, cz, eid});
    index_to_eid.push_back(eid);
    // Pre-allocate the bucket for this eid
    points_per_eid[eid] = {};
  }

  // Build the Index
  CentroidKDTree kdtree(3, centroids_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  kdtree.buildIndex();

  // Assign tree points to nearest edge
  size_t num_points = tree.size();
  for (size_t i = 0; i < num_points; ++i)
  {
    double query_pt[3] = { tree.get_x(i), tree.get_y(i), tree.get_z(i) };
    size_t ret_index;
    double out_dist_sqr;

    nanoflann::KNNResultSet<double> resultSet(1);
    resultSet.init(&ret_index, &out_dist_sqr);
    kdtree.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters());

    int eid = index_to_eid[ret_index];
    points_per_eid[eid].push_back(i);
  }

  return points_per_eid;
}

}

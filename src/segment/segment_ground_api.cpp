#include "arbor.h"
#include "PTD.h"
#include "myomp.h"
#include "nanoflann.h"

#include <vector>
#include <cmath>

namespace arbor::segment {

struct GroundCloud
{
  struct Point { double x, y, z; };
  std::vector<Point> pts;

  inline size_t kdtree_get_point_count() const { return pts.size(); }
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
    if (dim == 0) return pts[idx].x;
    return pts[idx].y;
  }
  template <class BBOX> bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};


void segment_ground(PointCloud& scene, const settings::ArborParameters& params, const Logger& logger)
{
  logger("Ground point segmentation");

  // Initialize bounds
  double xmin = std::numeric_limits<double>::max(), xmax = -std::numeric_limits<double>::max();
  double ymin = std::numeric_limits<double>::max(), ymax = -std::numeric_limits<double>::max();
  double zmin = std::numeric_limits<double>::max(), zmax = -std::numeric_limits<double>::max();
  for (unsigned int i = 0; i < scene.size() ; ++i)
  {
    double x = scene.get_x(i);
    double y = scene.get_y(i);
    double z = scene.get_z(i);

    if (x < xmin) xmin = x; if (x > xmax) xmax = x;
    if (y < ymin) ymin = y; if (y > ymax) ymax = y;
    if (z < zmin) zmin = z; if (z > zmax) zmax = z;
  }

  // Register the bounding box of the points
  PTD::Bbox bb;
  bb.xmin = xmin;
  bb.ymin = ymin;
  bb.zmin = zmin;
  bb.xmax = xmax;
  bb.ymax = ymax;
  bb.zmax = zmax;

  PTD::Parameters p;
  p.max_iteration_angle = 30;
  p.max_iteration_distance = 2;
  p.spacing = 0.1;
  p.seed_resolution_search = 5;
  p.verbose = true;

  PTD::PTD ptd(p, bb);

  for (unsigned int i = 0 ; i < scene.size() ; i++)
  {
    ptd.insert_point(scene.get_x(i), scene.get_y(i), scene.get_z(i), i);
  }

  ptd.run();

  std::vector<unsigned int> gnd_fids = ptd.get_ground_fid();
  GroundCloud ground_pts;
  ground_pts.pts.reserve(gnd_fids.size());

  for (size_t i = 0; i < scene.size() ; i++)
  {
    scene.is_ground(i);
    scene.set_classification(i, 1);
  }

  for (const auto fid : gnd_fids)
  {
    scene.set_classification(fid, 2);
    ground_pts.pts.push_back({scene.get_x(fid), scene.get_y(fid), scene.get_z(fid)});
  }

  logger("Height above ground");

  typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, GroundCloud>,
    GroundCloud, 2> my_kd_tree;

  my_kd_tree index(2, ground_pts, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  int nan_count = 0;

  #pragma omp parallel for reduction(+:nan_count)
  for (int i = 0; i < (int)scene.size(); ++i)
  {
    double cx = scene.get_x(i);
    double cy = scene.get_y(i);
    double cz = ptd.get_z(cx, cy);

    // Check if PTD returned NaN
    if (std::isnan(cz))
    {
      nan_count++;

      size_t ret_index;
      double out_dist_sqr;
      nanoflann::KNNResultSet<double> resultSet(1);
      resultSet.init(&ret_index, &out_dist_sqr);

      double query_pt[2] = {cx, cy};
      if (index.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParameters()))
      {
        cz = ground_pts.pts[ret_index].z;
      }
      else
      {
        cz = 0.0; // fallback
      }
    }
    scene.set_hag(i, scene.get_z(i) - cz);
  }

  if (nan_count > 0)
    logger("  " + std::to_string(nan_count) + " points outside TIN interpolated with 1-nn");
}

}

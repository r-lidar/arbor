#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

#include "nanoflann/nanoflann.h"

using namespace Rcpp;

// ============================================================================
//  1. Data Structures
// ============================================================================

struct Point3D
{
  double x, y, z;
};

// Simplified QSMcylinder based on your provided structure
struct QSMcylinder
{
  double startX, startY, startZ;
  double endX, endY, endZ;
  double radius;
  int cyl_ID;

  // Precomputed properties for optimization
  double ab_x, ab_y, ab_z; // Vector AB
  double ab_sqnorm;        // Squared length of AB
  double mid_x, mid_y, mid_z; // Midpoint for KD-tree

  void precompute()
  {
    ab_x = endX - startX;
    ab_y = endY - startY;
    ab_z = endZ - startZ;
    ab_sqnorm = ab_x*ab_x + ab_y*ab_y + ab_z*ab_z;
    mid_x = (startX + endX) / 2.0;
    mid_y = (startY + endY) / 2.0;
    mid_z = (startZ + endZ) / 2.0;
  }
};

// ============================================================================
//  2. Nanoflann Adaptor for Cylinder Midpoints
// ============================================================================

struct CylinderCloud
{
  std::vector<QSMcylinder> cylinders;

  inline size_t kdtree_get_point_count() const { return cylinders.size(); }

  inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
    if (dim == 0) return cylinders[idx].mid_x;
    if (dim == 1) return cylinders[idx].mid_y;
    return cylinders[idx].mid_z;
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<
  nanoflann::L2_Simple_Adaptor<double, CylinderCloud>,
  CylinderCloud, 3
> CylinderKDTree;

// ============================================================================
//  3. Geometric Logic
// ============================================================================

inline void get_point_cylinder_dist(const Point3D& p, const QSMcylinder& cyl, double& out_dist, double& out_r, int& out_id)
{
  // Vector AP = P - A
  double ap_x = p.x - cyl.startX;
  double ap_y = p.y - cyl.startY;
  double ap_z = p.z - cyl.startZ;

  // t = dot(AP, AB) / ||AB||^2
  double dot = ap_x * cyl.ab_x + ap_y * cyl.ab_y + ap_z * cyl.ab_z;
  double t = dot / cyl.ab_sqnorm;

  // Clamp t to segment [0, 1]
  if (t < 0.0) t = 0.0;
  else if (t > 1.0) t = 1.0;

  // Closest point C = A + t * AB
  double c_x = cyl.startX + t * cyl.ab_x;
  double c_y = cyl.startY + t * cyl.ab_y;
  double c_z = cyl.startZ + t * cyl.ab_z;

  // Euclidean distance ||P - C||
  double dx = p.x - c_x;
  double dy = p.y - c_y;
  double dz = p.z - c_z;
  double dist_center = std::sqrt(dx*dx + dy*dy + dz*dz);

  // Surface distance
  out_dist = std::max(0.0, dist_center - cyl.radius);
  out_r    = cyl.radius;
  out_id   = cyl.cyl_ID;
}

// ============================================================================
//  4. Main Rcpp Export
// ============================================================================

DataFrame compute_qsm_distances(DataFrame qsm_df, DataFrame pts_df)
{
  // --- 1. Parse QSM Data ---
  // Extract vectors from DataFrame
  NumericVector sX = qsm_df["startX"];
  NumericVector sY = qsm_df["startY"];
  NumericVector sZ = qsm_df["startZ"];
  NumericVector eX = qsm_df["endX"];
  NumericVector eY = qsm_df["endY"];
  NumericVector eZ = qsm_df["endZ"];
  NumericVector rad = qsm_df["radius"];
  IntegerVector id = qsm_df["cyl_ID"];

  int n_cyl = sX.size();
  CylinderCloud cloud;
  cloud.cylinders.reserve(n_cyl);

  for(int i = 0; i < n_cyl; ++i)
  {
    QSMcylinder c;
    c.startX = sX[i]; c.startY = sY[i]; c.startZ = sZ[i];
    c.endX = eX[i];   c.endY = eY[i];   c.endZ = eZ[i];
    c.radius = rad[i];
    c.cyl_ID = id[i];
    c.precompute(); // Calcs midpoints and vectors
    cloud.cylinders.push_back(c);
  }

  // --- 2. Build KD-Tree on Cylinders ---
  CylinderKDTree index(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  // --- 3. Parse Point Cloud Data ---
  NumericVector pX = pts_df["X"];
  NumericVector pY = pts_df["Y"];
  NumericVector pZ = pts_df["Z"];
  int n_pts = pX.size();

  // --- 4. Prepare Outputs ---
  NumericVector out_dist(n_pts);
  IntegerVector out_id(n_pts);
  NumericVector out_rad(n_pts);

  // --- 5. Computation Loop ---
  // Adaptive k: 1% of cylinders or 10, whichever is smaller (min 1)
  size_t k = std::min((size_t)10, (size_t)std::ceil(n_cyl * 0.01));
  if (k < 1) k = 1;

  std::vector<uint32_t> ret_index(k);
  std::vector<double> out_dist_sqr(k);

  for(int i = 0; i < n_pts; ++i)
  {
    Point3D p = { pX[i], pY[i], pZ[i] };
    double query_pt[3] = { p.x, p.y, p.z };

    // Find k nearest cylinder midpoints
    index.knnSearch(query_pt, k, &ret_index[0], &out_dist_sqr[0]);

    double min_d = std::numeric_limits<double>::max();
    int best_id = -1;
    double best_r = 0.0;

    // Check exact geometric distance for candidates
    for(size_t j = 0; j < k; ++j)
    {

      int idx = ret_index[j];
      double d_val, r_val;
      int id_val;

      get_point_cylinder_dist(p, cloud.cylinders[idx], d_val, r_val, id_val);

      if (d_val < min_d)
      {
        min_d = d_val;
        best_id = id_val;
        best_r = r_val;
      }
    }

    out_dist[i] = min_d;
    out_id[i]   = best_id;
    out_rad[i]  = best_r;
  }

  // Return as new DataFrame
  return DataFrame::create(
    Named("dist") = out_dist,
    Named("cyl_ID") = out_id,
    Named("radius") = out_rad
  );
}

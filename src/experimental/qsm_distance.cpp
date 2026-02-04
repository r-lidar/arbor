// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

#include "myomp.h"
#include "progressbar.h"
#include "nanoflann/nanoflann.h"

using namespace Rcpp;

// ============================================================================
//  1. Data Structures
// ============================================================================

struct QSMcylinder
{
  double startX, startY, startZ;
  double ab_x, ab_y, ab_z;   // Vector AB
  double ab_sqnorm;          // Squared length of AB
  double mid_x, mid_y, mid_z;
  double radius;
  int cyl_ID;

  void setup(double sX, double sY, double sZ, double eX, double eY, double eZ, double r, int id)
  {
    startX = sX; startY = sY; startZ = sZ;
    radius = r;
    cyl_ID = id;

    ab_x = eX - sX;
    ab_y = eY - sY;
    ab_z = eZ - sZ;
    ab_sqnorm = ab_x*ab_x + ab_y*ab_y + ab_z*ab_z;

    // Safety check for zero-length cylinders to avoid NaN
    if(ab_sqnorm < 1e-12) ab_sqnorm = 1e-12;

    mid_x = (sX + eX) * 0.5;
    mid_y = (sY + eY) * 0.5;
    mid_z = (sZ + eZ) * 0.5;
  }
};

// ============================================================================
//  2. Nanoflann Adaptor
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
//  3. Geometric Logic (Optimized)
// ============================================================================

// Returns true if a closer distance was found and updates best_*
inline void update_if_closer(const double px, const double py, const double pz,
                             const QSMcylinder& cyl,
                             double& best_dist, int& best_id, double& best_rad)
{
  // 1. Calculate projection t
  double ap_x = px - cyl.startX;
  double ap_y = py - cyl.startY;
  double ap_z = pz - cyl.startZ;

  double dot = ap_x * cyl.ab_x + ap_y * cyl.ab_y + ap_z * cyl.ab_z;
  double t = dot / cyl.ab_sqnorm;

  // 2. Clamp t to [0, 1] - manual logic is faster than std::clamp in some compilers
  if (t < 0.0) t = 0.0;
  else if (t > 1.0) t = 1.0;

  // 3. Distance to axis squared (avoid sqrt yet)
  double c_x = cyl.startX + t * cyl.ab_x;
  double c_y = cyl.startY + t * cyl.ab_y;
  double c_z = cyl.startZ + t * cyl.ab_z;

  double dx = px - c_x;
  double dy = py - c_y;
  double dz = pz - c_z;
  double dist_sq_axis = dx*dx + dy*dy + dz*dz;

  // 4. Optimization: "Lazy Sqrt"
  // We want to know if:  max(0, dist_axis - R) < best_dist
  // Logic:
  // If dist_axis < R (point inside cylinder), surface dist is 0.
  //    -> 0 is always <= best_dist (since best_dist >= 0). Update.
  // If dist_axis > R, we check: dist_axis - R < best_dist
  //    -> dist_axis < best_dist + R
  //    -> dist_axis^2 < (best_dist + R)^2

  double dist_check = best_dist + cyl.radius;

  // Only calculate sqrt if there is a chance this is the closest point
  if (dist_sq_axis < dist_check * dist_check) {
    double dist_axis = std::sqrt(dist_sq_axis);
    double surf_dist = (dist_axis > cyl.radius) ? (dist_axis - cyl.radius) : 0.0;

    if (surf_dist < best_dist) {
      best_dist = surf_dist;
      best_id   = cyl.cyl_ID;
      best_rad  = cyl.radius;
    }
  }
}

// ============================================================================
//  4. Main Rcpp Export
// ============================================================================

// [[Rcpp::export]]
DataFrame qsm_distances_cpp(DataFrame qsm_df, DataFrame pts_df)
{
  // --- 1. Parse QSM Data (Direct Pointer Access) ---
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
  cloud.cylinders.resize(n_cyl);

  // Standard loop for setup (this is fast enough usually)
  for(int i = 0; i < n_cyl; ++i) {
    cloud.cylinders[i].setup(sX[i], sY[i], sZ[i], eX[i], eY[i], eZ[i], rad[i], id[i]);
  }

  // --- 2. Build KD-Tree ---
  CylinderKDTree index(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  // --- 3. Parse Point Cloud Data ---
  NumericVector pX = pts_df["X"];
  NumericVector pY = pts_df["Y"];
  NumericVector pZ = pts_df["Z"];
  int n_pts = pX.size();

  // Prepare Output (std::vector for thread safety)
  std::vector<double> out_dist(n_pts);
  std::vector<int>    out_id(n_pts);
  std::vector<double> out_rad(n_pts);

  // Get raw pointers for thread-safe read access without Rcpp overhead
  const double* ptr_pX = &pX[0];
  const double* ptr_pY = &pY[0];
  const double* ptr_pZ = &pZ[0];

  size_t k = std::min((size_t)10, (size_t)std::ceil(n_cyl * 0.01));
  if (k < 1) k = 1;

  Progress pb(n_pts, "Point2Cylinders");
  std::atomic<bool> abort(false);

  #pragma omp parallel
  {
    // Thread-local storage to avoid reallocation inside loop
    std::vector<uint32_t> ret_index(k);
    std::vector<double> out_dist_sqr(k);

    #pragma omp for schedule(static)
    for(int i = 0; i < n_pts; ++i)
    {
      if (abort.load(std::memory_order_relaxed)) continue;
      if(pb.check_interrupt()) abort = true;
      pb.tick();

      double px = ptr_pX[i];
      double py = ptr_pY[i];
      double pz = ptr_pZ[i];
      double query_pt[3] = { px, py, pz };

      // 1. KNN Search
      index.knnSearch(query_pt, k, &ret_index[0], &out_dist_sqr[0]);

      // 2. Exact Distance Check
      double min_d = std::numeric_limits<double>::max();
      int best_id_val = -1;
      double best_r_val = 0.0;

      for(size_t j = 0; j < k; ++j)
      {
        int idx = ret_index[j];
        // Pass current min_d. The function will only update if it finds something smaller.
        update_if_closer(px, py, pz, cloud.cylinders[idx], min_d, best_id_val, best_r_val);
      }

      out_dist[i] = min_d;
      out_id[i]   = best_id_val;
      out_rad[i]  = best_r_val;
    }

    if (abort.load()) Rcpp::stop("Computation aborted");
  }

  // Return as new DataFrame
  return DataFrame::create(
    Named("dist") = wrap(out_dist),
    Named("cyl_ID") = wrap(out_id),
    Named("radius") = wrap(out_rad)
  );
}

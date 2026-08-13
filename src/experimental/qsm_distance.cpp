/**
 * @file qsm_distance.cpp
 * Project: Arbor
 *
 * Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

#include "myomp.h"
#include "services.h"
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

  inline double kdtree_get_pt(const size_t idx, const size_t dim) const
  {
    if (dim == 0) return cylinders[idx].mid_x;
    if (dim == 1) return cylinders[idx].mid_y;
    return cylinders[idx].mid_z;
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, CylinderCloud>, CylinderCloud, 3> CylinderKDTree;

// ============================================================================
//  3. Geometric Logic (Optimized)
// ============================================================================

// Returns true if a closer distance was found and updates best_*
inline void update_if_closer(const double px, const double py, const double pz,
                             const QSMcylinder& cyl,
                             double& best_dist, int& best_id, double& best_rad)
{
  // 1. Vector AP (Start to Point)
  double ap_x = px - cyl.startX;
  double ap_y = py - cyl.startY;
  double ap_z = pz - cyl.startZ;

  // 2. Project AP onto AB to find t
  double dot = ap_x * cyl.ab_x + ap_y * cyl.ab_y + ap_z * cyl.ab_z;
  double t = dot / cyl.ab_sqnorm;

  // Case A: The point projects onto the finite axis (0 <= t <= 1)
  // Distance is distance to the cylindrical surface
  if (t >= 0.0 && t <= 1.0)
  {
    double axis_pt_x = cyl.startX + t * cyl.ab_x;
    double axis_pt_y = cyl.startY + t * cyl.ab_y;
    double axis_pt_z = cyl.startZ + t * cyl.ab_z;

    double dx = px - axis_pt_x;
    double dy = py - axis_pt_y;
    double dz = pz - axis_pt_z;
    double dist_sq_to_axis = dx*dx + dy*dy + dz*dz;
    double dist_to_axis = std::sqrt(dist_sq_to_axis);

    // If inside the radius, distance is 0. If outside, it's (dist - R)
    double d = (dist_to_axis > cyl.radius) ? (dist_to_axis - cyl.radius) : 0.0;

    if (d < best_dist)
    {
      best_dist = d;
      best_id = cyl.cyl_ID;
      best_rad = cyl.radius;
    }
  }
  // Case B: The point projects beyond the ends (Caps)
  else
  {
    // Determine which cap we are closest to (Start or End)
    double cap_x, cap_y, cap_z;

    if (t < 0.0)
    {
      // Closest to Start Cap
      cap_x = cyl.startX;
      cap_y = cyl.startY;
      cap_z = cyl.startZ;
    }
    else
    {
      // t > 1.0
      // Closest to End Cap
      cap_x = cyl.startX + cyl.ab_x;
      cap_y = cyl.startY + cyl.ab_y;
      cap_z = cyl.startZ + cyl.ab_z;
    }

    // Vector from Cap Center to Point
    //double cp_x = px - cap_x;
    //double cp_y = py - cap_y;
    //double cp_z = pz - cap_z;

    // We need to find the distance to the *disk* defined by the cap.
    // However, since we are in the "beyond ends" region, the closest point
    // on the cylinder is either:
    // 1. On the rim of the cap (if point is outside the cylinder's infinite tube radius)
    // 2. On the flat face of the cap (if point is inside the infinite tube radius)

    // Calculate distance from point to the axis line (infinite)
    // We can reuse the projection logic or cross product.
    // Let's use the vector rejection from the axis.

    // Re-project onto axis relative to the specific cap is redundant if we trust t logic,
    // but let's be geometric:
    // Distance along the axis (longitudinal distance)
    double dist_along_axis;
    if (t < 0.0) dist_along_axis = std::sqrt(cyl.ab_sqnorm) * (-t); // Distance behind start
    else         dist_along_axis = std::sqrt(cyl.ab_sqnorm) * (t - 1.0); // Distance beyond end

    // Radial distance (distance from the infinite line)
    // We can compute the point on the infinite line corresponding to t
    double axis_inf_x = cyl.startX + t * cyl.ab_x;
    double axis_inf_y = cyl.startY + t * cyl.ab_y;
    double axis_inf_z = cyl.startZ + t * cyl.ab_z;

    double rad_dx = px - axis_inf_x;
    double rad_dy = py - axis_inf_y;
    double rad_dz = pz - axis_inf_z;
    double dist_radial = std::sqrt(rad_dx*rad_dx + rad_dy*rad_dy + rad_dz*rad_dz);

    double final_dist;

    if (dist_radial <= cyl.radius)
    {
      // Point is within the tube's radius, but past the end.
      // Distance is purely the longitudinal distance to the flat cap face.
      final_dist = dist_along_axis;
    }
    else
    {
      // Point is outside the tube and past the end.
      // The closest point is on the circular rim (edge of the cap).
      // We form a right triangle:
      // leg 1: distance from cap plane (dist_along_axis)
      // leg 2: distance from rim in radial direction (dist_radial - radius)
      double d_rim = dist_radial - cyl.radius;
      final_dist = std::sqrt(dist_along_axis*dist_along_axis + d_rim*d_rim);
    }

    if (final_dist < best_dist)
    {
      best_dist = final_dist;
      best_id = cyl.cyl_ID;
      best_rad = cyl.radius;
    }
  }
}

// ============================================================================
//  4. Main Rcpp
// ============================================================================

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

  for(int i = 0; i < n_cyl; ++i)
  {
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

  // Prepare Output
  std::vector<double> out_dist(n_pts);
  std::vector<int>    out_id(n_pts);
  std::vector<double> out_rad(n_pts);

  // Get raw pointers for thread-safe read access without Rcpp overhead
  const double* ptr_pX = &pX[0];
  const double* ptr_pY = &pY[0];
  const double* ptr_pZ = &pZ[0];

  size_t k = std::min((size_t)10, (size_t)std::ceil(n_cyl * 0.01));
  if (k < 1) k = 1;

  auto pb = ServiceLocator::make_progress(n_pts, "Point2Cylinders");
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
      if(pb->check_interrupt()) abort = true;
      pb->tick();

      double px = ptr_pX[i];
      double py = ptr_pY[i];
      double pz = ptr_pZ[i];
      double query_pt[3] = { px, py, pz };

      // 1. KNN Search
      index.knnSearch(query_pt, k, &ret_index[0], &out_dist_sqr[0]);

      // 2. Exact Distance Check
      double min_d = std::numeric_limits<double>::max();
      int best_id_val = NA_INTEGER;
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

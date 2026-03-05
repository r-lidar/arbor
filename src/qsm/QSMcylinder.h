#ifndef QSMCYLINDER_H
#define QSMCYLINDER_H

#include <cmath>

namespace arbor::qsm {

static constexpr double SUBTREE_LENGTH_UNSET    = -1.0;
static constexpr double SUBTREE_MAXZ_UNSET      = -1e300;
static constexpr double SUBTREE_VOLUME_UNSET    = -1;
static constexpr double RADIUS_UNSET            = -1.0;
static constexpr double Z_EPS                   = 1e-9;

struct QSMcylinder
{
  double startX = 0;
  double startY = 0;
  double startZ = 0;
  double endX = 0;
  double endY = 0;
  double endZ = 0;
  double radius = RADIUS_UNSET;
  double conic_allometry = RADIUS_UNSET;
  double subtree_length = SUBTREE_LENGTH_UNSET;
  double subtree_max_endZ = SUBTREE_MAXZ_UNSET;
  double subtree_volume = SUBTREE_VOLUME_UNSET;
  int cyl_ID = 0;
  int parent_ID = 0;
  int axis_ID = 0;
  int branch_order = 0;

  inline double length() const noexcept
  {
    double dx = endX - startX;
    double dy = endY - startY;
    double dz = endZ - startZ;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  inline double volume() const noexcept
  {
    double r2 = radius * radius;
    return M_PI * r2 * length();
  }

  // Compute distance on the fly
  inline double distance(double px, double py, double pz) const
  {
    // Compute Vector AB (Cylinder Axis)
    double ab_x = endX - startX;
    double ab_y = endY - startY;
    double ab_z = endZ - startZ;

    // Compute Squared Length of AB
    double ab_sq = ab_x*ab_x + ab_y*ab_y + ab_z*ab_z;
    if(ab_sq < 1e-12) ab_sq = 1e-12; // Avoid division by zero

    // Compute Vector AP (Start to Point)
    double ap_x = px - startX;
    double ap_y = py - startY;
    double ap_z = pz - startZ;

    // Project AP onto AB (Dot Product)
    double t = (ap_x*ab_x + ap_y*ab_y + ap_z*ab_z) / ab_sq;

    // Clamp t to [0, 1] to stay within the segment
    if (t < 0.0) t = 0.0;
    else if (t > 1.0) t = 1.0;

    // Find Closest Point on Axis
    double c_x = startX + t * ab_x;
    double c_y = startY + t * ab_y;
    double c_z = startZ + t * ab_z;

    // Distance from Point to Axis
    double dx = px - c_x;
    double dy = py - c_y;
    double dz = pz - c_z;

    double dist_axis = std::sqrt(dx*dx + dy*dy + dz*dz);

    // Subtract radius for surface distance (clamped to 0)
    return (dist_axis > radius) ? (dist_axis - radius) : 0.0;
  }

  inline double angle()
  {
    double dx = endX - startX;
    double dy = endY - startY;
    double dz = endZ - startZ;
    double len = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (len < 1e-6) return 0.0;
    return std::acos(dz / len) * 180.0 / M_PI;
  }
};

}

#endif

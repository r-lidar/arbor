#include "QSM.h"

#include <map>
#include <cmath>
#include <algorithm>

#include "QSM.h"
#include "nanoflann/Adaptor.h"
#include "nanoflann/nanoflann.h"
#include "ransac.h"

using PointCloud = DataFrameAdaptor;
typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, SimpleAdaptor>, SimpleAdaptor, 3> CentroidKDTree;

inline double QSM::conic_allometry(double tip_radius, double wi, double w0, double r0) const
{
  const double s = std::pow(wi / w0, 1.1);
  return tip_radius + s * (r0 - tip_radius);
}

/********************************/

void QSM::measure_radii(const PointCloud& tree)
{
  if (cylinders_.empty())
    throw std::runtime_error("Internal error: no cylinder in this QSM");

  // Prepare centroids for KD-Tree
  SimpleAdaptor centroids_cloud;
  centroids_cloud.points.reserve(cylinders_.size());

  // We also need a way to map the KD-tree index back to the cylinder pointer efficiently
  std::vector<QSMcylinder*> index_to_cyl_ptr;
  index_to_cyl_ptr.reserve(cylinders_.size());

  for (auto& pair : cylinders_)
  {
    QSMcylinder& cyl = pair.second;
    double cx = (cyl.startX + cyl.endX) / 2.0;
    double cy = (cyl.startY + cyl.endY) / 2.0;
    double cz = (cyl.startZ + cyl.endZ) / 2.0;
    centroids_cloud.points.push_back({cx, cy, cz, cyl.cyl_ID});
    index_to_cyl_ptr.push_back(&cyl);
  }

  // Build KD-Tree on Centroids
  CentroidKDTree index(3, centroids_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  // Assign tree points to nearest cylinder (Edge of the QSM graph)
  // Structure: Vector of vectors, where outer index matches 'index_to_cyl_ptr'
  std::vector<std::vector<size_t>> points_per_cylinder(centroids_cloud.points.size());
  size_t num_points = tree.point_count();

  for (size_t i = 0; i < num_points; ++i)
  {
    double query_pt[3];
    query_pt[0] = tree.get_x(i);
    query_pt[1] = tree.get_y(i);
    query_pt[2] = tree.get_z(i);

    size_t ret_index;
    double out_dist_sqr;

    nanoflann::KNNResultSet<double> resultSet(1);
    resultSet.init(&ret_index, &out_dist_sqr);
    index.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters());

    points_per_cylinder[ret_index].push_back(i);
  }

  // Iterate over cylinders and run RANSAC
  for (size_t i = 0; i < points_per_cylinder.size(); ++i)
  {
    const auto& point_indices = points_per_cylinder[i];
    QSMcylinder* cyl = index_to_cyl_ptr[i];

    if (cyl->cyl_ID < 1) continue;              // cyl_ID < 1 means it is a prolongation.
    if (point_indices.empty()) continue;        // No point for this cylinder
    if (cyl->conic_allometry < 0.04) continue;  // Skip if the conic allometry predicts a too small radius
    if (point_indices.size() < 50) continue;    // Skip ransac if we have very few points.

    // Compute ransdac circle for this cylinder
    // providing the orientation of the cylinder for alignment
    RansacCircle ransac(100, 0.02);
    for (size_t pt_idx : point_indices)
    {
      double x = tree.get_x(pt_idx);
      double y = tree.get_y(pt_idx);
      double z = tree.get_z(pt_idx);
      ransac.add_point(x,y,z);
    }
    std::array<double, 3> axis_start = {cyl->startX, cyl->startY, cyl->startZ};
    std::array<double, 3> axis_end   = {cyl->endX,   cyl->endY,   cyl->endZ};
    ransac.find_circle(axis_start, axis_end);

    // Validate we have good circle fitting
    double r_meas = ransac.get_radius();
    double p_inside = ransac.get_inside_percentage(); // 0.0 to 1.0
    double arc = ransac.get_arc_coverage();           // degrees
    double p_inlier = ransac.get_inlier_percentage(); // 0.0 to 1.0
    bool valid = true;
    if (r_meas < 0.05) valid = false;
    else if (p_inside > 0.20) valid = false;
    else if (!(arc > 180.0 && p_inlier > 0.30)) valid = false;

    // If the circle is valid then we have a measurement
    // for this cylinder
    if (valid) cyl->radius = r_meas;
  }
}

void QSM::polynomial_fitting(double tip_radius)
{
  // Group Cylinder IDs by Axis ID
  std::map<int, std::vector<int>> axes;
  for (auto& kv : cylinders_)
  {
    const QSMcylinder& cyl = kv.second;
    axes[cyl.axis_ID].push_back(cyl.cyl_ID);
  }

  // Iterate over each axis
  for (auto& group : axes)
  {
    int axis_id = group.first;
    const std::vector<int>& cyl_ids = group.second;

    // Collect valid measurements for regression
    // We need to solve: (radius - tip) = a * length + b * length^2
    // Let y = radius - tip_radius
    // Let x = subtree_length
    // Equation: y = a*x + b*x^2

    std::vector<std::pair<double, double>> data_points;
    data_points.reserve(cyl_ids.size());

    for (int id : cyl_ids)
    {
      // Access via map to ensure we get the actual object
      const auto& cyl = cylinders_[id];

      // Corresponds to R: sum(!is.na(axe$radius))
      if (cyl.radius != RADIUS_UNSET && cyl.subtree_length != SUBTREE_LENGTH_UNSET)
      {
        double y = cyl.radius - tip_radius;
        double x = cyl.subtree_length;
        data_points.push_back({x, y});
      }
    }

    // No interpolation if less than 6 points
    if (data_points.size() <= 6)
    {
      continue;
    }

    // 3. Solve Linear Least Squares (OLS) for 2 variables (a, b)
    // We want to minimize sum of squared errors.
    // System of Normal Equations for y = c1*x + c2*x^2:
    // | sum(x^2)  sum(x^3) |  | a |   | sum(x*y)   |
    // | sum(x^3)  sum(x^4) |  | b | = | sum(x^2*y) |

    double sum_x2 = 0.0;
    double sum_x3 = 0.0;
    double sum_x4 = 0.0;
    double sum_xy = 0.0;
    double sum_x2y = 0.0;

    for (const auto& p : data_points)
    {
      double x = p.first;
      double y = p.second;
      double x2 = x * x;
      double x3 = x2 * x;
      double x4 = x3 * x;

      sum_x2  += x2;
      sum_x3  += x3;
      sum_x4  += x4;
      sum_xy  += x * y;
      sum_x2y += x2 * y;
    }

    // Determinant of the 2x2 matrix
    double det = sum_x2 * sum_x4 - sum_x3 * sum_x3;

    // Check for singular matrix (collinear points or not enough spread)
    // Fallback or skip: Matrix is singular, cannot fit polynomial.
    if (std::abs(det) < 1e-12)
    {
      continue;
    }

    // Solve using Cramer's rule
    double a = (sum_xy * sum_x4 - sum_x2y * sum_x3) / det;
    double b = (sum_x2 * sum_x2y - sum_xy * sum_x3) / det;


    // Apply prediction to ALL cylinders in this axis
    for (int id : cyl_ids)
    {
      QSMcylinder& cyl = cylinders_[id];

      if (cyl.subtree_length != SUBTREE_LENGTH_UNSET)
      {
        double len = cyl.subtree_length;
        double new_radius = tip_radius + a * len + b * len * len;
        if (new_radius < 0) new_radius = 0; // Safety clamp: Radius shouldn't be negative
        cyl.radius = new_radius;
      }
    }
  }
}

void QSM::reconstruct_missing_radii(double tip_radius)
{
  // Group cylinders by "branch_order" for efficient iteration.
  std::map<int, std::vector<QSMcylinder*>> cylinders_by_branch_order;
  for (auto& [_, cyl] : cylinders_)  cylinders_by_branch_order[cyl.branch_order].push_back(&cyl);

  // Loop by branch order. Start at 2, main trunk (order 1) should already have been computed
  for (auto& [branch_order, cyls] : cylinders_by_branch_order)
  {
    if (branch_order == 1) continue;

    // Group cylinders of the current "branch_order" by "axis_ID"
    std::unordered_map<int, std::vector<QSMcylinder*>> axes;
    for (QSMcylinder* c : cyls) axes[c->axis_ID].push_back(c);

    // Loop on each axis
    for (auto& [axis_id, axe] : axes)
    {
      if (axe.empty()) continue;

      // Check if ANY cylinder in the axis needs reconstruction. This means that polynomial
      // fitting failed.
      // Iterate through cylinders to check for UNSET radius and count valid measurements
      bool needs_reconstruction = false;
      for (const QSMcylinder* c : axe)
      {
        if (c->radius == RADIUS_UNSET)
        {
          needs_reconstruction = true;
          break;
        }
      }
      if (!needs_reconstruction) continue;

      // Sort the cylinder by ID
      std::sort(axe.begin(), axe.end(), [](const QSMcylinder* a, const QSMcylinder* b) { return a->cyl_ID < b->cyl_ID; });

      // The first cylinder is the root of the branch. We search for its parent
      const int parent_id = axe.front()->parent_ID;
      QSMcylinder& parent = get_cylinder_by_id(parent_id);
      const double r0 = parent.radius;
      const double w0 = parent.subtree_length;

      // Compute theoretical radii by conic radiometry. This is our fallback value
      std::vector<double> r_theoretical;
      r_theoretical.reserve(axe.size());
      for (const QSMcylinder* c : axe)
        r_theoretical.push_back(conic_allometry(tip_radius, c->subtree_length, w0, r0));

      // Collect measurements if any
      std::vector<double> r_measured(axe.size());
      int valid_measurements = 0;

      std::transform(axe.begin(), axe.end(), r_measured.begin(), [&](const QSMcylinder* c)
      {
        if (c->radius != RADIUS_UNSET)
        {
          valid_measurements++;
          return c->radius;
        }
        return std::numeric_limits<double>::quiet_NaN(); // preserve alignment
      });

      // By default the final measurement are the theoretical values
      std::vector<double> r_final = r_theoretical;

      // If we have more than 3 measurements we can try to estimate a better profile
      if (valid_measurements >= 3)
      {
        // Compute ratios only for valid measured radii
        std::vector<double> ratios;
        ratios.reserve(valid_measurements);
        for (size_t k = 0; k < axe.size(); ++k)
        {
          if (!std::isnan(r_measured[k]))
            ratios.push_back(r_measured[k] / r_theoretical[k]);
        }

        // Median ratio calculation
        std::sort(ratios.begin(), ratios.end());
        const size_t m = ratios.size() / 2;
        const double ratio_median =  (ratios.size() % 2 != 0) ? ratios[m]: (ratios[m - 1] + ratios[m]) / 2.0;

        // Apply scaling using the conic allometry formulation
        for (size_t k = 0; k < axe.size(); ++k)
        {
          const double wi = axe[k]->subtree_length;
          r_final[k] = conic_allometry(tip_radius, wi, w0, r0 * ratio_median);
        }
      }

      // Update QSM
      for (size_t k = 0; k < axe.size(); ++k)
      {
        axe[k]->radius = r_final[k];
      }
    }
  }
}


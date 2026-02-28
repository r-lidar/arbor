#include "QSM.h"

#include <map>
#include <cmath>
#include <algorithm>

#include "QSM.h"
#include "PointCloud.h"
#include "nanoflann.h"
#include "ransac.h"

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

inline double QSM::conic_allometry(double tip_radius, double wi, double w0, double r0) const
{
  const double s = std::pow(wi / w0, 1.1);
  return tip_radius + s * (r0 - tip_radius);
}

/********************************/

// @param sarc sensitivity arc
// @param sins sensitiviy inside
// @param sinl sensitiviy inliner
// @param srmeas sensitivity r measured
void QSM::measure_radii(const PointCloud& tree, float sarc, float sins, float sinl, float srmeas)
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
    if (r_meas < srmeas) valid = false;
    else if (p_inside > sins) valid = false;
    else if (!(arc > sarc && p_inlier > sinl)) valid = false;

    // If the circle is valid then we have a measurement
    // for this cylinder
    if (valid) cyl->radius = r_meas;
  }
}

void QSM::polynomial_fitting(double tip_radius)
{
  // Group cylindera by axis
  std::map<int, Axe> axes;
  for (auto& kv : cylinders_)
  {
    QSMcylinder& cyl = kv.second;
    axes[cyl.axis_ID].add_cylinder(&cyl);
  }

  // Iterate over each axis
  for (auto& group : axes)
  {
    int axis_id = group.first;
    Axe& axe = group.second;
    axe.sort();

    // Collect valid measurements for regression
    // We need to solve: (radius - tip) = a * length + b * length^2
    // Let y = radius - tip_radius
    // Let x = subtree_length
    // Equation: y = a*x + b*x^2

    std::vector<std::pair<double, double>> data_points;
    data_points.reserve(axe.size());

    for (QSMcylinder* cyl : axe)
    {
      // Corresponds to R: sum(!is.na(axe$radius))
      if (cyl->radius != RADIUS_UNSET && cyl->subtree_length != SUBTREE_LENGTH_UNSET)
      {
        double y = cyl->radius - tip_radius;
        double x = cyl->subtree_length;
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
    if (std::abs(det) < 1e-12) continue;

    // Solve using Cramer's rule
    double a = (sum_xy * sum_x4 - sum_x2y * sum_x3) / det;
    double b = (sum_x2 * sum_x2y - sum_xy * sum_x3) / det;

    // Compute the whole prediction for this axis
    std::vector<double> predictions;
    for (QSMcylinder* cyl : axe)
    {
      double len = cyl->subtree_length;
      double pred_radius = tip_radius + a * len + b * len * len;
      if (pred_radius < 0) pred_radius = 0; // Safety clamp
      predictions.push_back(pred_radius);
    }

    // Sometime the polynomial fitting can mess-up the branch with parent cylinder
    // smaller than child. Usually this happens if the measure we have are centered
    // on the branch with no data on the insertion point (or root). In this case the
    // polynomial is not constrained and weird behavior can arise. So we check that
    // we have no increasing diameters and repair otherwise.
    double previous_radius = predictions.back();
    for (auto it = predictions.rbegin(); it != predictions.rend(); ++it)
    {
      if (*it < previous_radius) { *it = previous_radius; }
      previous_radius = *it;
    }

    // Apply prediction to cylinders in this axis
    int i = 0;
    for (QSMcylinder* cyl : axe)
    {
      double pred_radius = predictions[i]; i++;

      if (cyl->subtree_length == SUBTREE_LENGTH_UNSET)
        throw std::logic_error("subtree_length unset during polynomial fitting");

      bool should_update = false;

      // Condition 1: Value is missing
      if (cyl->radius == RADIUS_UNSET)
      {
        should_update = true;
      }
      // Condition 2: Measured radius is > X% different than prediction
      else
      {
        // Calculate absolute difference
        double diff = 99999; // Desactivate feature
        //double diff = std::abs(cyl->radius - pred_radius);

        // Check if difference > X% of the PREDICTED value
        // (Using prediction as the baseline for the "ideal" curve)
        if (diff > 0.20 * pred_radius)
        {
          should_update = true;
        }
      }

      if (should_update)
      {
        cyl->radius = pred_radius;
      }
    }
  }
}

void QSM::reconstruct_missing_radii(double tip_radius)
{
  // Group cylinders by "branch_order"
  std::map<int, std::vector<QSMcylinder*>> cylinders_by_branch_order;
  for (auto& [_, cyl] : cylinders_)  cylinders_by_branch_order[cyl.branch_order].push_back(&cyl);

  // Loop by branch order. Start at 2, main trunk (order 1) should already have been computed
  for (auto& [branch_order, cyls] : cylinders_by_branch_order)
  {
    if (branch_order == 1) continue;

    // Group cylinders of the current "branch_order" by "axis_ID". Order does not matter
    std::unordered_map<int, Axe> axes;
    for (QSMcylinder* c : cyls) axes[c->axis_ID].add_cylinder(c);

    // Loop on each axis
    for (auto& [axis_id, axe] : axes)
    {
      if (axe.empty()) continue;

      // Sort the cylinder by subtree_length to get them from root to tip
      axe.sort();

      // The first cylinder is the root of the branch. We search for its parent
      const int parent_id = axe[0]->parent_ID;
      QSMcylinder& parent = get_cylinder_by_id(parent_id);
      const double r0 = parent.radius*0.9;
      const double w0 = parent.subtree_length;

      // If the axes has missing radii it means no polynomial fitting was performed
      // we need to reconstruct the radii
      if (axe.need_reconstruction())
      {
        // Compute theoretical radii by conic allometry. This is our fallback value
        for (QSMcylinder* c : axe)
          c->radius = conic_allometry(tip_radius, c->subtree_length, w0, r0);
      }
      // If the axes has no missing radii it means polynomial fitting was performed
      // However sometime the fitting may generate branch bigger than their parent. This
      // fixes such ugly output.
      else
      {
        double r1 = axe[0]->radius;
        double ratio = r0/r1;
        if (ratio < 1)  axe.scale(ratio);
      }
    }
  }
}


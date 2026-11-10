/**
 * @file segment_ground_api.cpp
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

#include "arbor.h"
#include "PTD.h"
#include "myomp.h"
#include "nanoflann.h"
#include "progress.h"

#include <algorithm>
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


// Quantile of the values (linear interpolation between the closest ranks, like R's default quantile type 7).
// The values are reordered. They must not be empty.
static float quantile(std::vector<float>& values, double probability)
{
  const size_t n = values.size();

  // Fractional rank of the quantile, e.g. 2.4 means 40% of the way between the 3rd and 4th smallest values
  const double rank       = (n - 1) * probability;
  const size_t lower_rank = static_cast<size_t>(std::floor(rank));
  const double fraction   = rank - lower_rank;

  // Partial sort: the value at lower_rank is in place, all the values after it are larger or equal
  auto lower_it = values.begin() + lower_rank;
  std::nth_element(values.begin(), lower_it, values.end());
  const float lower_value = *lower_it;

  if (lower_rank + 1 >= n) return lower_value;

  // The next rank is the smallest value after lower_it
  const float upper_value = *std::min_element(lower_it + 1, values.end());
  return lower_value + fraction * (upper_value - lower_value);
}

// Median of the values. The values are reordered. They must not be empty.
static float median(std::vector<float>& values)
{
  return quantile(values, 0.5);
}

// Pre-classifies as noise (7) the points far below the ground, before the ground segmentation, so they cannot be
// picked as ground. Low noise is common in TLS/MLS (multipath, reflections) and is the lowest thing in the scene:
// a ground filter would happily take it as the ground.
//
// 1. A grid of resolution noise_res stores, in each cell, a low quantile of Z (noise_quantile): a rough ground.
// 2. This rough ground is smoothed with a median over noise_window x noise_window cells: the reference ground.
// 3. The points lower than noise_threshold below the reference ground are classified noise.
//
// Only the eligible points are used. Only the eligible points that are not already classified (e.g. ground set by
// the user) can become noise; those are then no longer eligible. Returns the number of points classified noise.
static size_t classify_low_noise(PointCloud& scene, std::vector<bool>& eligible, const settings::GroundParameters& p)
{
  const size_t n = scene.size();
  const int NOISE_CLASS = 7;

  // ---- Grid covering the eligible points ----

  double xmin =  std::numeric_limits<double>::max();
  double xmax = -std::numeric_limits<double>::max();
  double ymin =  std::numeric_limits<double>::max();
  double ymax = -std::numeric_limits<double>::max();

  for (size_t i = 0 ; i < n ; ++i)
  {
    if (!eligible[i]) continue;
    const double x = scene.get_x(i);
    const double y = scene.get_y(i);
    xmin = std::min(xmin, x);
    xmax = std::max(xmax, x);
    ymin = std::min(ymin, y);
    ymax = std::max(ymax, y);
  }

  const bool no_eligible_point = xmin > xmax;
  if (no_eligible_point) return 0;

  // A raster must have a non-zero extent
  if (xmax == xmin) xmax += p.noise_res;
  if (ymax == ymin) ymax += p.noise_res;

  // ---- Step 1: rough ground, a low quantile of Z in each cell ----

  Raster rough_ground(xmin, ymin, xmax, ymax, p.noise_res);
  const size_t ncells = rough_ground.get_ncells();

  // Z of the eligible points, grouped by cell
  std::vector<std::vector<float>> z_by_cell(ncells);
  for (size_t i = 0 ; i < n ; ++i)
  {
    if (!eligible[i]) continue;
    const size_t cell = rough_ground.cell_from_xy(scene.get_x(i), scene.get_y(i));
    z_by_cell[cell].push_back(scene.get_z(i));
  }

  const double probability = std::clamp(p.noise_quantile, 0.0, 1.0);

  #pragma omp parallel for
  for (size_t cell = 0 ; cell < ncells ; ++cell)
  {
    std::vector<float>& z = z_by_cell[cell];
    if (z.empty()) continue; // the cell stays NA

    rough_ground.set_value(cell, quantile(z, probability));
    std::vector<float>().swap(z); // free the memory as soon as possible
  }

  // ---- Step 2: reference ground, the rough ground smoothed with a focal median ----

  const Raster reference_ground = rough_ground.focal(p.noise_window, median);

  // ---- Step 3: points far below the reference ground are noise ----

  size_t n_noise = 0;
  for (size_t i = 0 ; i < n ; ++i)
  {
    const bool can_be_reclassified = eligible[i] && scene.is_processable(i);
    if (!can_be_reclassified) continue;

    const float ground_z = reference_ground.get_value(scene.get_x(i), scene.get_y(i));
    if (std::isnan(ground_z)) continue;

    const double height = scene.get_z(i) - ground_z;
    const bool is_below_ground = height < -p.noise_threshold;
    if (!is_below_ground) continue;

    scene.set_classification(i, NOISE_CLASS);
    eligible[i] = false;
    n_noise++;
  }

  return n_noise;
}

void segment_ground(PointCloud& scene, const settings::ArborParameters& params)
{
  if (!scene.has_class()) throw std::runtime_error("segment_ground: classification atttibute not allocated.");

  size_t n = scene.size();
  if (n == 0) throw std::runtime_error("segment_ground: point cloud is empty.");

  Task task("Ground segmentation", {1, 8, 2});

  // Initialize bounds
  double xmin = std::numeric_limits<double>::max(), xmax = -std::numeric_limits<double>::max();
  double ymin = std::numeric_limits<double>::max(), ymax = -std::numeric_limits<double>::max();
  double zmin = std::numeric_limits<double>::max(), zmax = -std::numeric_limits<double>::max();

  // Use pre-classified data to exclude points. Works only if we read pre-classified data.
  // Only the candidates for the segmentation and the points already classified ground are classified.
  // The other points (noise, buildings...) keep their class: a low noise point must never become ground.
  std::vector<bool> eligible(n);
  size_t n_eligible = 0;
  for (size_t i = 0; i < n ; ++i)
  {
    eligible[i] = scene.is_processable(i) || (scene.has_class() && scene.get_classification(i) == 2);
    if (eligible[i]) n_eligible++;
  }
  if (n_eligible == 0) throw std::runtime_error("segment_ground: no point to classify (only noise, buildings...).");

  task.step("Low noise classification");

  // The general case is that the previous code block does nothing because the point cloud is not classified
  // We pre-classify below ground noise points
  if (params.ground.noise_enabled)
  {
    size_t n_noise = classify_low_noise(scene, eligible, params.ground);
    if (n_noise > 0) logging::info(std::to_string(n_noise) + " points below the ground classified noise");
    n_eligible -= n_noise;
    if (n_eligible == 0) throw std::runtime_error("segment_ground: no point to classify (only noise, buildings...).");
  }

  task.step("Ground point classification");

  for (unsigned int i = 0; i < scene.size() ; ++i)
  {
    if (!eligible[i]) continue;
    double x = scene.get_x(i);
    double y = scene.get_y(i);
    double z = scene.get_z(i);

    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
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

  for (unsigned int i = 0 ; i < n ; i++)
  {
    if (!eligible[i]) continue;
    ptd.insert_point(scene.get_x(i), scene.get_y(i), scene.get_z(i), i);
  }

  ptd.run();

  std::vector<unsigned int> gnd_fids = ptd.get_ground_fid();
  GroundCloud ground_pts;
  ground_pts.pts.reserve(gnd_fids.size());

  for (size_t i = 0; i < n ; i++)
  {
    if (eligible[i]) scene.set_classification(i, 1);
  }

  for (const auto fid : gnd_fids)
  {
    scene.set_classification(fid, 2);
    ground_pts.pts.push_back({scene.get_x(fid), scene.get_y(fid), scene.get_z(fid)});
  }

  task.step("Height above ground");

  // HAG is computed from the TIN of ground points (ptd.get_z). The TIN is only defined inside the
  // convex hull of the ground points, so points outside (typically at the edges of the scene) get NaN.
  // For those we fall back to the z of the nearest ground point in 2D. The KD-tree is built on the
  // ground points only and is queried only for these NaN cases.
  
  typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, GroundCloud>, GroundCloud, 2> my_kd_tree;

  my_kd_tree index(2, ground_pts, nanoflann::KDTreeSingleIndexAdaptorParams(10, nanoflann::KDTreeSingleIndexAdaptorFlags::SkipInitialBuildIndex));
  index.buildIndex();

  int nan_count = 0;

  #pragma omp parallel for reduction(+:nan_count)
  for (size_t i = 0; i < n; ++i)
  {
    double cx = scene.get_x(i);
    double cy = scene.get_y(i);
    double cz = ptd.get_z(cx, cy);

    // NaN means (cx, cy) is outside the TIN: use the nearest ground point (1-NN) instead
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
    double z = scene.get_z(i);
    scene.set_hag(i, z - cz);
  }

  if (nan_count > 0)
    logging::info(std::to_string(nan_count) + " points outside TIN interpolated with 1-nn");
}

}

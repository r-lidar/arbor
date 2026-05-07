/**
 * @file SeedDetector.cpp
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

#include "SeedDetector.h"
#include "Grid3D.h"

namespace arbor::seeds {

void SeedDetector::run(const PointCloud& scene)
{
  auto log = ServiceLocator::logger();

  if (scene.size() == 0)     throw std::runtime_error("find_seed: point cloud is empty.");
  if (!scene.has_hag())      throw std::runtime_error("find_seed: point cloud is missing required 'hag' attribute.");
  if (!scene.has_foliage())  throw std::runtime_error("find_seed: point cloud is missing required 'foliage' attribute.");

  log("Seed detection start");

  log("Computing lower bound");
  compute_min_hag(scene);

  log("Extracting passages");
  extract_passages(scene);

  log("Slicing the point cloud");
  slice_wood(scene);

  log("Circle detection");
  circles = detect_tree_circles(wood);
  if (circles.empty()) throw std::runtime_error("No circle detected in wood slices");

  log("Generate tree cages");
  make_cages();

  log("Safe zone exclusion");
  safe_zone();

  log("Find primary seeds");
  find_primary_seeds();

  log("Pathfinder");
  merge_short_passages();

  log("Find secondary seeds");
  find_secondary_seeds();

  seeds = primary_seeds + secondary_seeds;

  log("Seed detection completed");
}

void SeedDetector::compute_min_hag(const PointCloud& scene)
{
  min_hag = scene.get_hag(0);
  for (size_t i = 1 ; i < scene.size() ; i++)
  {
    if (scene.get_hag(i) < min_hag)
      min_hag = scene.get_hag(i);
  }
}

void SeedDetector::extract_passages(const PointCloud& scene)
{
  size_t n = scene.size();

  // The heights at which we slice
  // We extract some slices of wood (thickness ~3cm) at different
  // height including min_hag + input parameters typically 25, 70, 90 cm
  double thick = params.seeds.slice_thickness;
  std::vector<double> heights = params.seeds.slice_at;
  heights.push_back(min_hag+thick);
  std::sort(heights.begin(), heights.end());
  double th = heights.back() + 0.2;

  // Long passage selection below a threshold height
  std::vector<bool> long_passage_mask(n, false);
  for (size_t i = 0 ; i < n ; i++)
    long_passage_mask[i] = (scene.get_passage(i) > params.seeds.min_passage) && (scene.get_hag(i) < th);

  // Short passage selection 50 cm above the cut
  std::vector<bool> short_passage_mask(n, false);
  for (size_t i = 0 ; i < n ; i++)
    short_passage_mask[i] = (scene.get_passage(i) > 1) && (scene.get_passage(i) < 10) && (scene.get_hag(i) < min_hag + 0.5);

  long_passages  = scene.subset(long_passage_mask);
  short_passages = scene.subset(short_passage_mask);

  long_passages  = densify_passages(long_passages);
  short_passages = densify_passages(short_passages);
}

PointCloud SeedDetector::densify_passages(const PointCloud& x)
{
  PointCloud xup = x;
  PointCloud xdw = x;
  xup.translate(0, 0, 0.01);
  xdw.translate(0, 0, -0.01);
  for (size_t i = 0 ; i < x.size() ; i++)
  {
    xup.set_passage(i, -1);
    xdw.set_passage(i, -1);
  }

  PointCloud y = x + xup + xdw;
  return y;
}

// Extract slices of wood points below a threshold
// -----------------------------------------------
void SeedDetector::slice_wood(const PointCloud& scene)
{
  size_t n = scene.size();

  // The heights at which we slice
  // We extract some slices of wood (thickness ~3cm) at different
  // height including min_hag + input parameters typically 25, 70, 90 cm
  double thick = params.seeds.slice_thickness;
  std::vector<double> heights = params.seeds.slice_at;
  heights.push_back(min_hag+thick);
  std::sort(heights.begin(), heights.end());

  std::vector<bool> wood_mask(n, false);
  double half = params.seeds.slice_thickness * 0.5;

  for (std::size_t i = 0; i < n; ++i)
  {
    if (!scene.is_wood(i)) continue;

    double h = scene.get_hag(i);

    for (const double s : heights)
    {
      if (h > (s - half) && h < (s + half))
      {
        wood_mask[i] = true;
        break;
      }
    }
  }

  wood = scene.subset(wood_mask);
}

// Safe zone exclusion
// We exclude some wood points where we have circles
// This allows to separate contact trees where we detected
// circles
// -------------------------------------------------
void SeedDetector::safe_zone()
{
  const double safe_zone = 0.2;
  const double inner_threshold = 0.02;

  std::vector<bool> keep(wood.size(), true);

  // For each point in the wood point cloud
  for (size_t i = 0; i < wood.size(); i++)
  {
    double px = wood.get_x(i);
    double py = wood.get_y(i);
    double pz = wood.get_z(i);

    // Check distance to each circle
    for (size_t j = 0; j < circles.size(); j++)
    {
      double cx = circles[j].X;
      double cy = circles[j].Y;
      double cz = circles[j].Z;
      double r = circles[j].R;

      // Calculate 3D distance from point to circle center
      double d = std::sqrt((px - cx) * (px - cx) + (py - cy) * (py - cy) + (pz - cz) * (pz - cz));

      // Mark point for removal if it's in the safe zone
      if (d > (r + inner_threshold) && d < (r + safe_zone))
      {
        keep[i] = false;
        break; // No need to check other circles for this point
      }
    }
  }

  // Create subset of wood without the safe zone points
  wood = wood.subset(keep);
}

// Main tree seeds
// Merging the cage, the wood and the passage we compute
// connected components
// -------------------------------------------------
void SeedDetector::find_primary_seeds()
{
  // Merge wood slices + long passage of the path finder + cages
  PointCloud temp = wood + long_passages + cages;

  // Compute connected components of that point cloud to assign an ID to each point
  temp.scale(1, 1, 0.5);
  double res = std::round(params.pathfinder.decimation * 100.0) / 100.0;
  Grid3D grid(temp, res);
  std::vector<int> cluster_ids = grid.connected_components(26);
  for (size_t i = 0 ; i < temp.size() ; i++) temp.set_treeid(i, cluster_ids[i]);
  temp.scale(1, 1, 2);

  // We retain only the point with a passage number > 0
  // This effectively removes wood points that are subset upstream without passages (assigned 0)
  // It retains cages are assigned passage = 9999
  std::vector<bool> long_seed_mask(temp.size(), false);
  for (size_t i = 0; i < temp.size(); i++) long_seed_mask[i] = (temp.get_passage(i) > 0);

  primary_seeds = temp.subset(long_seed_mask);
}

void SeedDetector::merge_short_passages()
{
  if (short_passages.size() == 0) return;
  
  // Force short passage to be wood to avoid wood/foliage penalties in pathfinder
  for (size_t i = 0; i < short_passages.size(); i++) short_passages.set_foliage(i, 0);

  // Setup parameters for pathfinder
  settings::ArborParameters p;
  p.pathfinder.max_gap = 0.1;
  p.pathfinder.k = 10;
  p.pathfinder.k_seed = 2;
  p.pathfinder.power = 1;
  p.pathfinder.angle_penalty = std::vector<float>(181);
  std::fill(p.pathfinder.angle_penalty.begin(), p.pathfinder.angle_penalty.end(), 1.0f);

  arbor::segment::segment_instance(short_passages, primary_seeds, p);

  std::vector<bool> has_id_mask(short_passages.size(), false);
  for (size_t i = 0; i < short_passages.size(); i++)
  {
    int tree_id = short_passages.get_treeid(i);
    has_id_mask[i] = (tree_id >= 0); // -1 or negative means NA
  }

  primary_seeds += short_passages.subset(has_id_mask);
}

void SeedDetector::find_secondary_seeds()
{
  std::vector<bool> no_id_mask(short_passages.size(), false);
  for (size_t i = 0; i < short_passages.size(); i++)
  {
    int tree_id = short_passages.get_treeid(i);
    no_id_mask[i] = (tree_id < 0); // Assuming -1 or negative means NA
  }

  PointCloud short_passages_noid = short_passages.subset(no_id_mask);

  int max_id = primary_seeds.get_treeid(0);
  for (size_t i = 1 ; i < primary_seeds.size() ; i++)
  {
    if (primary_seeds.get_treeid(i) > max_id)
      max_id = primary_seeds.get_treeid(i);
  }

  if (short_passages_noid.size() > 0)
  {
    Grid3D grid(short_passages_noid, 0.1);
    std::vector<int> short_cluster_ids = grid.connected_components(26);

    for (size_t i = 0; i < short_passages_noid.size(); i++)
      short_passages_noid.set_treeid(i, short_cluster_ids[i] + max_id + 1);
  }

  secondary_seeds = short_passages_noid;
}

void SeedDetector::filter_seeds()
{
  // Retain only seeds below breast height (1m)
  std::vector<bool> mask(seeds.size(), false);
  for (size_t i = 0; i < seeds.size(); i++) mask[i] = (seeds.get_hag(i) < 1 && seeds.get_passage(i) > 0);
  seeds = seeds.subset(mask);
}

}

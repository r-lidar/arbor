/**
 * @file qsf_api.cpp
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
#include "myomp.h"

#include <map>
#include <vector>
#include <string>
#include <atomic>

namespace arbor::qsm {

QSF qsf(const PointCloud& scene, double min_height, const settings::ArborParameters& params)
{
  if (min_height < 2)      throw std::runtime_error("Height limit cannot be inferior to 2");
  if (!scene.has_hag())    throw std::runtime_error("Missing attribute 'hag' in the point cloud");
  if (!scene.has_treeid()) throw std::runtime_error("Missing attribute 'treeID' in the point cloud");
  if (scene.size() == 0)   throw std::runtime_error("Point cloud with 0 point: failure");

  QSF result;

  // Group indices and track max height per ID
  std::unordered_map<int, std::vector<unsigned int>> tree_indices;
  std::unordered_map<int, double> tree_heights;

  for (size_t i = 0; i < scene.size(); ++i)
  {
    int id = scene.get_treeid(i);
    if (id <= 0 || id == std::numeric_limits<int>::max()) continue; // NA from R or negative (removed) or 0 from arbor

    int userdata = scene.get_userdata(i);
    if (userdata != 0) continue; // Not ARBORTREE

    double hag = scene.get_hag(i);

    tree_indices[id].push_back(i);

    if (tree_heights.find(id) == tree_heights.end() || hag > tree_heights[id])
    {
      tree_heights[id] = hag;
    }
  }

  // Filter IDs that meet the height requirement immediately
  std::vector<int> valid_tree_ids;
  valid_tree_ids.reserve(tree_indices.size());
  for (auto const& [id, height] : tree_heights)
  {
    if (height >= min_height)
      valid_tree_ids.push_back(id);
  }
  std::sort(valid_tree_ids.begin(), valid_tree_ids.end());

  auto p = ServiceLocator::make_progress(valid_tree_ids.size(), "QSF");
  auto old_logger = ServiceLocator::logger();
  ServiceLocator::register_logger([](const std::string&) {});
  std::atomic<bool> error_occurred{false};
  std::exception_ptr eptr = nullptr;

  #pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < valid_tree_ids.size(); ++i)
  {
    // Check if another thread encountered an error
    if (error_occurred.load(std::memory_order_relaxed)) continue;

    int current_id = valid_tree_ids[i];
    const std::vector<unsigned int>& indices = tree_indices[current_id];

    try
    {
      PointCloud tree = scene.subset(indices);
      QSM q = qsm(tree, params);
      q.id = current_id;
      q.name = "tree_" + std::to_string(current_id);

      #pragma omp critical
      {
        result.add_qsm(q);
      }
    }
    catch (...)
    {
      #pragma omp critical
      {
        if (!error_occurred)
        {
          error_occurred = true;
          try
          {
            std::rethrow_exception(std::current_exception());
          }
          catch (const std::exception& e)
          {
            eptr = std::make_exception_ptr(
              std::runtime_error("Tree ID " + std::to_string(current_id) + ": " + e.what())
            );
          }
          catch (...)
          {
            eptr = std::make_exception_ptr(
              std::runtime_error("Tree ID " + std::to_string(current_id) + ": unknown error")
            );
          }
        }
      }
    }

    p->tick();
  }

  // If something went wrong, rethrow the exception outside the parallel region
  if (eptr)
  {
    ServiceLocator::register_logger(old_logger);
    std::rethrow_exception(eptr);
  }

  p->finalize();
  ServiceLocator::register_logger(old_logger);
  return result;
}

}

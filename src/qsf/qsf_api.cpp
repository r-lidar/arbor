#include "arbor.h"
#include "QSF.h"
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
  std::unordered_map<int, std::vector<int>> tree_indices;
  std::unordered_map<int, double> tree_heights;

  for (size_t i = 0; i < scene.size(); ++i)
  {
    int id = scene.get_treeid(i);
    if (id < 0) continue; // NA from R or -1 from arbor
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
    const std::vector<int>& indices = tree_indices[current_id];

    try
    {
      PointCloud tree = scene.subset(indices);
      QSM q = qsm(tree, params);

      #pragma omp critical
      {
        result.add_qsm(std::to_string(current_id), q);
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

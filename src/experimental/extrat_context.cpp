#include "nanoflann.h"
#include "Adaptor.h"

#include <unordered_set>
#include <vector>

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 3> KDTree;

std::vector<int> extract_tree_context(const PointCloud& pc, int target_tree_id,  bool exclude_tree = false, int k = 10)
{
  // Find all points belonging to the target tree
  std::vector<size_t> target_indices;
  target_indices.reserve(pc.size() / 100);
  for (size_t i = 0; i < pc.size(); ++i)
  {
    if (pc.get_treeid(i) == target_tree_id)
    {
      target_indices.push_back(i);
    }
  }

  if (target_indices.empty())
  {
    throw std::runtime_error("Requested tree ID is not part of the point cloud");
  }

  KDTree index(3, pc, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  std::unordered_set<int> contact_tree_ids;

  std::vector<unsigned int> ret_index(k);
  std::vector<double> out_dist_sqr(k);
  for (size_t target_idx : target_indices)
  {
    double query_pt[3] = {pc.get_x(target_idx), pc.get_y(target_idx), pc.get_z(target_idx)};

    size_t n_found = index.knnSearch(query_pt, k, &ret_index[0], &out_dist_sqr[0]);

    for (size_t j = 0; j < n_found; ++j)
    {
      size_t neighbor_idx = ret_index[j];
      int neighbor_tree_id = pc.get_treeid(neighbor_idx);

      if (neighbor_tree_id > 0 && neighbor_tree_id != target_tree_id)
      {
        contact_tree_ids.insert(neighbor_tree_id);
      }
    }
  }

  std::vector<int> result;
  if (!exclude_tree)
  {
    result.push_back(target_tree_id);
  }

  result.insert(result.end(), contact_tree_ids.begin(), contact_tree_ids.end());
  std::sort(result.begin(), result.end());

  return result;
}

// Rcpp wrapper
Rcpp::IntegerVector extract_tree_context_cpp(Rcpp::DataFrame las, int tree_id,  bool exclude_tree = false, int k = 10)
{
  PointCloud adaptor(las);
  std::vector<int> result = extract_tree_context(adaptor, tree_id, exclude_tree, k);
  return Rcpp::wrap(result);
}

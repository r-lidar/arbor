#ifndef POINTCLOUDBASE_H
#define POINTCLOUDBASE_H

#include <cmath>
#include <memory>
#include <vector>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <random>
#include <algorithm>

#include "nanoflann.h"

inline void hcl_to_rgb(float h, float c, float l, uint8_t* R, uint8_t* G, uint8_t* B)
{
  // 1. Convert HCL to CIELAB
  float h_rad = h * M_PI / 180.0f;
  float L = l;
  float a = std::cos(h_rad) * c;
  float b = std::sin(h_rad) * c;

  // 2. Convert CIELAB to XYZ
  auto f_inv = [](float t)
  {
    return (t > 6.0f/29.0f) ? (t * t * t) : (3.0f * (6.0f/29.0f) * (6.0f/29.0f) * (t - 4.0f/29.0f));
  };

  float y = (L + 16.0f) / 116.0f;
  float x = y + a / 500.0f;
  float z = y - b / 200.0f;

  // Scale by D65 white point
  x = 0.95047f * f_inv(x);
  y = 1.00000f * f_inv(y);
  z = 1.08883f * f_inv(z);

  // 3. Convert XYZ to Linear RGB
  float r_lin =  3.2406f * x - 1.5372f * y - 0.4986f * z;
  float g_lin = -0.9689f * x + 1.8758f * y + 0.0415f * z;
  float b_lin =  0.0557f * x - 0.2040f * y + 1.0570f * z;

  // 4. Gamma correction (sRGB) and Clamping
  auto gamma = [](float val)
  {
    val = std::max(0.0f, std::min(1.0f, val));
    return (val <= 0.0031308f) ? (12.92f * val) : (1.055f * std::pow(val, 1.0f/2.4f) - 0.055f);
  };

  *R = gamma(r_lin)*255;
  *G = gamma(g_lin)*255;
  *B = gamma(b_lin)*255;
}

// ============================================================================
// DECLARATIONS (Interface)
// ============================================================================

// CRTP design pattern. BasePointCloud provides a common interface and implements
// functions that can operate on the derived class through compile-time polymorphism.
// The base class does not know the memory layout of the derived class, which is
// defined by the inherited classes.
//
// PointCloudDataFrame maps directly onto the memory of an R data.frame to integrate
// with an R package. PointCloudDefault manages its own memory layout and is used
// by ArborStudio.
template <typename Derived>
class BasePointCloud
{
  using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, Derived>, Derived, 3>;

protected:
  size_t n_points = 0;
  size_t true_n_points = 0;
  mutable std::unique_ptr<KDTree> kdtree;

  // 0. Default constructor
  BasePointCloud() = default;
  virtual ~BasePointCloud() = default;

  // 1. Custom Copy Constructor
  BasePointCloud(const BasePointCloud& other) : n_points(other.n_points), true_n_points(other.true_n_points), kdtree(nullptr) {}

  // 2. Custom Copy Assignment Operator
  BasePointCloud& operator=(const BasePointCloud& other) {
    if (this != &other) {
      n_points = other.n_points;
      true_n_points = other.true_n_points;
      kdtree = nullptr;
    }
    return *this;
  }

  // 3. Default Move Semantics (Transfers the unique_ptr safely)
  BasePointCloud(BasePointCloud&&) noexcept = default;
  BasePointCloud& operator=(BasePointCloud&&) noexcept = default;

  // Needed for partition
  virtual void swap_points(size_t i, size_t j) = 0;

public:
  enum class UserData : uint8_t
  {
    TREE       = 0,
    LOW        = 1,
    UNDERSTORY = 2,
    BUFFER     = 3
  };

  // For nanoflann
  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const;
  inline size_t kdtree_get_point_count() const;

  // Generic
  inline size_t size() const;
  inline size_t true_size() const;

  void build_index() const;
  void knn(const double* query, int k, std::vector<unsigned int>& indices, std::vector<double>& sqdist) const;
  void radius_search(const double* query, double radius, std::vector<unsigned int>& indices, std::vector<double>& sqdist) const;
  Derived subset(const std::vector<bool>& mask, bool xyz_only = false) const;
  Derived subset_by_treeid(int tid) const;
  std::vector<unsigned int> find_points_by_tree_id(int tid) const;
  std::vector<unsigned int> find_context_points_by_tree_id(int tid, float r = 0.1, bool exclude_tree = false) const;


  // Partition like std::partition
  template <typename Predicate>
  size_t partition(Predicate pred);

  void colorize_trees(bool darken_foliage = false);
};

template <typename Derived>
template <class BBOX>
bool BasePointCloud<Derived>::kdtree_get_bbox(BBOX&) const
{
  return false;
}

template <typename Derived>
inline size_t BasePointCloud<Derived>::kdtree_get_point_count() const
{
  return size();
}

template <typename Derived>
inline size_t BasePointCloud<Derived>::size() const
{
  return n_points;
}

template <typename Derived>
inline size_t BasePointCloud<Derived>::true_size() const
{
  return true_n_points;
}

template <typename Derived>
void BasePointCloud<Derived>::build_index() const
{
  if (kdtree) return;
  const Derived& self = static_cast<const Derived&>(*this);
  kdtree = std::make_unique<KDTree>(3, self, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  kdtree->buildIndex();
}

template <typename Derived>
void BasePointCloud<Derived>::knn(const double* query, int k, std::vector<unsigned int>& indices, std::vector<double>& sqdist) const
{
  if (!kdtree) build_index();

  indices.clear();
  sqdist.clear();

  if (k <= 0) return;

  indices.resize(k);
  sqdist.resize(k);

  nanoflann::KNNResultSet<double, unsigned int> resultSet(k);
  resultSet.init(indices.data(), sqdist.data());

  nanoflann::SearchParameters params;
  params.eps = 0.0f;
  params.sorted = false;

  kdtree->findNeighbors(resultSet, query, params);
}

template <typename Derived>
void BasePointCloud<Derived>::radius_search(const double* query, double radius, std::vector<unsigned int>& indices, std::vector<double>& sqdist) const
{
  if (!kdtree) build_index();

  indices.clear();
  sqdist.clear();

  if (radius <= 0.0) return;

  double radius_sq = radius * radius;

  std::vector<nanoflann::ResultItem<unsigned int, double>> ret_matches;

  nanoflann::SearchParameters params;
  params.eps = 0.0f;
  params.sorted = false;

  // Note: nanoflann natively expects the squared search radius
  const size_t nMatches = kdtree->radiusSearch(query, radius_sq, ret_matches, params);

  indices.reserve(nMatches);
  sqdist.reserve(nMatches);

  for (size_t i = 0; i < nMatches; ++i)
  {
    indices.push_back(ret_matches[i].first);
    sqdist.push_back(ret_matches[i].second);
  }
}

template <typename Derived>
Derived BasePointCloud<Derived>::subset(const std::vector<bool>& mask, bool xyz_only) const
{
  const Derived& self = static_cast<const Derived&>(*this);
  if (mask.size() != self.size())
    throw std::runtime_error("subset mask size mismatch: expected " + std::to_string(self.size()) + " but got " + std::to_string(mask.size()));

  std::vector<unsigned int> indices;
  indices.reserve(mask.size() / 10); // 10 % heuristic is arbitrary
  for (size_t i = 0; i < mask.size(); ++i)
  {
    if (mask[i]) indices.push_back(static_cast<unsigned int>(i));
  }

  if (indices.size() == self.size()) return self;
  return self.subset(indices, xyz_only);
}

template <typename Derived>
Derived BasePointCloud<Derived>::subset_by_treeid(int tid) const
{
  const Derived& self = static_cast<const Derived&>(*this);

  std::vector<int> keep;
  for (size_t i = 0 ; i < self.size() ; i++)
  {
    if (self.get_treeid(i) == tid)
      keep.push_back(i);
  }

  return self.subset(keep, false);
}

template <typename Derived>
std::vector<unsigned int> BasePointCloud<Derived>::find_points_by_tree_id(int tid) const
{
  const Derived& self = static_cast<const Derived&>(*this);

  std::vector<unsigned int> keep;
  for (size_t i = 0 ; i < self.size() ; i++)
  {
    if (self.get_treeid(i) == tid)
      keep.push_back(i);
  }

  return keep;
}

template <typename Derived>
template <typename Predicate>
size_t BasePointCloud<Derived>::partition(Predicate pred)
{
  if (n_points == 0) return 0;

  int64_t left = 0;
  int64_t right = static_cast<int64_t>(n_points) - 1;

  while (left <= right)
  {
    // Move left pointer until we find a point that should be at the back
    while (left <= right && pred(left))
    {
      left++;
    }
    // Move right pointer until we find a point that should be at the front
    while (left <= right && !pred(right))
    {
      right--;
    }

    if (left <= right)
    {
      swap_points(left, right);
      left++;
      right--;
    }
  }

  // Update n_points so nanoflann only sees the front section
  size_t new_size = static_cast<size_t>(left);
  n_points = new_size;

  kdtree.reset();

  return new_size;
}

template <typename Derived>
void BasePointCloud<Derived>::colorize_trees(bool darken_foliage)
{
  Derived& self = static_cast<Derived&>(*this);

  if (!self.has_red() || !self.has_green() || !self.has_blue())
    throw std::runtime_error("RGB memory not allocated");

  if (!self.has_treeid())
    throw std::runtime_error("No treeID in this point cloud");

  if (!self.has_foliage())
    darken_foliage = false;

  if (self.true_size() == 0) return;

  struct RGB { uint8_t r, g, b; };
  std::unordered_map<int, RGB> color_cache;
  const float darken_factor = 0.7f;

  for (size_t i = 0; i < self.size(); ++i)
  {
    int id = self.get_treeid(i);
    if (id <= 0) continue;

    auto [it, inserted] = color_cache.try_emplace(id, RGB{});
    if (inserted)
    {
      std::mt19937 gen(static_cast<uint32_t>(id));
      std::uniform_real_distribution<float> dist_h(0.0f, 360.0f);
      std::uniform_real_distribution<float> dist_c(42.0f, 98.0f);
      std::uniform_real_distribution<float> dist_l(40.0f, 90.0f);
      hcl_to_rgb(dist_h(gen), dist_c(gen), dist_l(gen),                 &it->second.r, &it->second.g, &it->second.b);
    }

    RGB color = it->second;

    // Abstracted wood/foliage check
    if (darken_foliage && !self.is_wood(i))
    {
      color.r *= darken_factor;
      color.g *= darken_factor;
      color.b *= darken_factor;
    }

    // Unified assignment
    self.set_color(i, color.r, color.g, color.b);
  }
}

template <typename Derived>
std::vector<unsigned int> BasePointCloud<Derived>::find_context_points_by_tree_id(int tid, float r, bool exclude_tree) const
{
  const Derived& self = static_cast<const Derived&>(*this);

  // Find all points belonging to the target tree
  std::vector<unsigned int> target_indices = self.find_points_by_tree_id(tid);

  if (target_indices.empty())
  {
    throw std::runtime_error("Requested tree ID is not part of the point cloud");
  }

  self.build_index();

  std::unordered_set<int> contact_tree_ids;

  std::vector<unsigned int> ret_index;
  std::vector<double> out_dist_sqr;
  for (size_t target_idx : target_indices)
  {
    double q[3];
    self.get_point(target_idx, q);

    self.radius_search(q, r, ret_index, out_dist_sqr);

    for (size_t j = 0; j < ret_index.size(); ++j)
    {
      size_t neighbor_idx = ret_index[j];
      int neighbor_tree_id = self.get_treeid(neighbor_idx);

      if (neighbor_tree_id > 0 && neighbor_tree_id < std::numeric_limits<int>::max() && neighbor_tree_id != tid)
      {
        contact_tree_ids.insert(neighbor_tree_id);
      }
    }
  }

  if (!exclude_tree)
  {
    contact_tree_ids.insert(tid);
  }

  std::vector<unsigned int> result;

  for (size_t i = 0 ; i < self.size() ; i++)
  {
    int tid = self.get_treeid(i);
    if (contact_tree_ids.count(tid) > 0)
      result.push_back(static_cast<unsigned int>(i));
  }

  std::sort(result.begin(), result.end());

  return result;
}


#endif

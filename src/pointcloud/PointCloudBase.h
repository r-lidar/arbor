#ifndef POINTCLOUDBASE_H
#define POINTCLOUDBASE_H

#include <cmath>
#include <memory>
#include <vector>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <algorithm>
#include <limits>

#include "nanoflann.h"
#include "Polygon.h"

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
  void remap_treeid(int from, int to);
  Derived subset(const std::vector<bool>& mask, bool xyz_only = false) const;
  Derived subset_by_treeid(int tid) const;

  std::vector<unsigned int> find_points_by_tree_id(int tid) const;
  std::vector<unsigned int> find_context_points_by_tree_id(int tid, float r = 0.1, bool exclude_tree = false) const;


  // Partition like std::partition
  template <typename Predicate>
  size_t partition(Predicate pred);

  void colorize_trees(bool darken_foliage = false);


  void flag_small_trees(double max_height);
  void flag_buffer(double buffer, const Derived* seeds = nullptr);
  void flag_buffer(const arbor::utils::Polygon& polygon, const Derived* seeds = nullptr);

private:
  struct Root { int treeid; double x, y; };
  std::vector<Root> compute_roots(const Derived* seeds) const;
};

// Same algorithm as R's mean() for doubles (long double accumulation + correction pass) so results
// are bit-identical to the former R implementation
template <typename Get>
inline double r_mean(size_t n, Get get)
{
  long double s = 0.0;
  for (size_t i = 0; i < n; ++i) s += get(i);
  s /= n;
  if (std::isfinite(static_cast<double>(s)))
  {
    long double t = 0.0;
    for (size_t i = 0; i < n; ++i) t += (get(i) - s);
    s += t / n;
  }
  return static_cast<double>(s);
}

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
  // Assigned once built: buildIndex() throws arbor::Cancelled on cancellation
  auto tree = std::make_unique<KDTree>(3, self, nanoflann::KDTreeSingleIndexAdaptorParams(10, nanoflann::KDTreeSingleIndexAdaptorFlags::SkipInitialBuildIndex));
  tree->buildIndex();
  kdtree = std::move(tree);
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
void BasePointCloud<Derived>::remap_treeid(int from, int to)
{
  Derived& self = static_cast<Derived&>(*this);

  if (from == to) return;

  if (!self.has_treeid())
    throw std::runtime_error("No treeID in this point cloud");

  for (size_t i = 0; i < self.size(); ++i)
  {
    if (self.get_treeid(i) == from)
    {
      self.set_treeid(i, to);
    }
  }
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

template <typename Derived>
void BasePointCloud<Derived>::flag_small_trees(double max_height)
{
  Derived& self = static_cast<Derived&>(*this);

  if (!self.has_treeid()) throw std::runtime_error("No treeID in this point cloud");
  if (!self.has_hag()) throw std::runtime_error("No hag in this point cloud");
  if (!self.has_userdata()) throw std::runtime_error("No UserData in this point cloud");

  // Height range of each tree. NaN hag are ignored
  struct Range { double min, max; };
  std::unordered_map<int, Range> ranges;
  constexpr double inf = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < self.size(); ++i)
  {
    int id = self.get_treeid(i);
    if (id <= 0) continue;

    double h = self.get_hag(i);
    if (std::isnan(h)) continue;
    auto [it, inserted] = ranges.try_emplace(id, Range{inf, -inf});
    Range& r = it->second;
    r.min = std::min(r.min, h);
    r.max = std::max(r.max, h);
  }

  std::unordered_set<int> keep;
  for (const auto& [id, r] : ranges)
  {
    if (r.max > max_height && r.min < max_height) keep.insert(id);
  }

  constexpr int NA = std::numeric_limits<int>::min();
  const int TREE = static_cast<int>(UserData::TREE);
  const int LOW = static_cast<int>(UserData::LOW);
  const int UNDERSTORY = static_cast<int>(UserData::UNDERSTORY);
  const int BUFFER = static_cast<int>(UserData::BUFFER);

  // BUFFER takes precedence over UNDERSTORY and is never overwritten
  for (size_t i = 0; i < self.size(); ++i)
  {
    int ud = self.get_userdata(i);
    if (ud == NA || ud == BUFFER) continue;
    if (ud == UNDERSTORY) ud = TREE;
    if (ud != LOW && keep.count(self.get_treeid(i)) == 0) ud = UNDERSTORY;
    self.set_userdata(i, ud);
  }
}

template <typename Derived>
std::vector<typename BasePointCloud<Derived>::Root> BasePointCloud<Derived>::compute_roots(const Derived* seeds) const
{
  const Derived& self = static_cast<const Derived&>(*this);
  std::vector<Root> roots;

  if (seeds)
  {
    if (!seeds->has_treeid()) throw std::runtime_error("No treeID in the seeds");

    // Mean position of the seeds of each tree
    std::unordered_map<int, std::vector<unsigned int>> groups;
    for (size_t i = 0; i < seeds->size(); ++i)
    {
      int id = seeds->get_treeid(i);
      if (id > 0) groups[id].push_back(static_cast<unsigned int>(i));
    }

    for (const auto& [id, idx] : groups)
    {
      double x = r_mean(idx.size(), [&](size_t k) { return seeds->get_x(idx[k]); });
      double y = r_mean(idx.size(), [&](size_t k) { return seeds->get_y(idx[k]); });
      roots.push_back({id, x, y});
    }
  }
  else
  {
    // Mean position of the 100 lowest points of each tree. Points are ordered by Z then by index
    // (stable sort) and NaN are last (same as R's order())
    constexpr size_t nmax = 100;
    constexpr size_t nmin = 10;

    // Z is stored with the index to avoid a cache miss when comparing to the top of the heap
    struct Item { double z; unsigned int i; };
    auto before = [](const Item& a, const Item& b)
    {
      bool na = std::isnan(a.z);
      bool nb = std::isnan(b.z);
      if (na != nb) return nb;
      if (!na && a.z != b.z) return a.z < b.z;
      return a.i < b.i;
    };

    // Max-heap of the nmax lowest points of each tree
    struct Group { size_t n = 0; std::vector<Item> heap; };
    std::unordered_map<int, Group> groups;
    for (size_t i = 0; i < self.size(); ++i)
    {
      int id = self.get_treeid(i);
      if (id <= 0) continue;

      Group& g = groups[id];
      g.n++;
      Item item{self.get_z(i), static_cast<unsigned int>(i)};
      if (g.heap.size() < nmax)
      {
        g.heap.push_back(item);
        std::push_heap(g.heap.begin(), g.heap.end(), before);
      }
      else if (before(item, g.heap.front()))
      {
        std::pop_heap(g.heap.begin(), g.heap.end(), before);
        g.heap.back() = item;
        std::push_heap(g.heap.begin(), g.heap.end(), before);
      }
    }

    for (auto& [id, g] : groups)
    {
      if (g.n < nmin) continue;
      std::sort_heap(g.heap.begin(), g.heap.end(), before);
      const std::vector<Item>& h = g.heap;
      double x = r_mean(h.size(), [&](size_t k) { return self.get_x(h[k].i); });
      double y = r_mean(h.size(), [&](size_t k) { return self.get_y(h[k].i); });
      roots.push_back({id, x, y});
    }
  }

  roots.erase(std::remove_if(roots.begin(), roots.end(), [](const Root& r) { return !std::isfinite(r.x) || !std::isfinite(r.y); }), roots.end());
  return roots;
}

template <typename Derived>
void BasePointCloud<Derived>::flag_buffer(const arbor::utils::Polygon& polygon, const Derived* seeds)
{
  Derived& self = static_cast<Derived&>(*this);

  if (!self.has_treeid()) throw std::runtime_error("No treeID in this point cloud");
  if (!self.has_userdata()) throw std::runtime_error("No UserData in this point cloud");

  std::unordered_set<int> rooted;
  std::unordered_set<int> keep;
  for (const Root& r : compute_roots(seeds))
  {
    rooted.insert(r.treeid);
    if (polygon.point_in_polygon(r.x, r.y)) keep.insert(r.treeid);
  }

  constexpr int NA = std::numeric_limits<int>::min();
  const int TREE = static_cast<int>(UserData::TREE);
  const int LOW = static_cast<int>(UserData::LOW);
  const int BUFFER = static_cast<int>(UserData::BUFFER);

  // LOW points are never flagged. Trees are flagged according to their root. Unsegmented points and
  // points of trees without a root (no seed, or too few points) are flagged according to their own
  // position. BUFFER overrides UNDERSTORY
  for (size_t i = 0; i < self.size(); ++i)
  {
    int ud = self.get_userdata(i);
    if (ud == NA || ud == LOW) continue;
    if (ud == BUFFER) ud = TREE;
    int id = self.get_treeid(i);
    bool inside = rooted.count(id) ? keep.count(id) > 0 : polygon.point_in_polygon(self.get_x(i), self.get_y(i));
    if (!inside) ud = BUFFER;
    self.set_userdata(i, ud);
  }
}

template <typename Derived>
void BasePointCloud<Derived>::flag_buffer(double buffer, const Derived* seeds)
{
  const Derived& self = static_cast<const Derived&>(*this);

  if (buffer > 0) throw std::invalid_argument("'buffer' must be negative");

  arbor::utils::Polygon hull = arbor::utils::Polygon::convex_hull(self.size(), [&](size_t i) { return self.get_x(i); }, [&](size_t i) { return self.get_y(i); });
  flag_buffer(hull.shrink(-buffer), seeds);
}

#endif

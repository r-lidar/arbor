#ifndef USING_R

#include "PointCloudDefault.h"

PointCloudDefault::PointCloudDefault()
{
}

PointCloudDefault::PointCloudDefault(size_t n, Attr attrs)
{
  try
  {
    coords.resize(n, {0.0f, 0.0f, 0.0f});
    n_points = n;
    allocate(attrs);
  }
  catch (...)
  {
    cleanup();
    throw;
  }
}

// ------------------------------------------------------------
// Insertion
// ------------------------------------------------------------
void PointCloudDefault::add_point(float x, float y, float z)
{
  if (n_points != coords.size())
    throw std::logic_error("add_point() cannot be called inside a ScopedPartition");

  coords.push_back({x,y,z});
  visit_attributes([this](Attr flag, auto member, auto def)
  {
    if (has(flag)) (this->*member).push_back(def);
  });
  n_points++;
  kdtree.reset();
}

// ------------------------------------------------------------
// Attribute management
// ------------------------------------------------------------
void PointCloudDefault::allocate(Attr a)
{
  // coords.size() rather than n_points: inside a ScopedPartition the hidden points need the attribute too
  const size_t n = coords.size();
  visit_attributes([this, a, n](Attr flag, auto member, auto def)
  {
    if (!any(a & flag) || has(flag)) return;
    (this->*member).assign(n, def);
    present |= flag;
  });
}

void PointCloudDefault::release(Attr a)
{
  visit_attributes([this, a](Attr flag, auto member, auto)
  {
    if (!any(a & flag)) return;
    auto& v = this->*member;
    v.clear();
    v.shrink_to_fit();
    present &= ~flag;
  });
}

// ------------------------------------------------------------
// Cleanup
// ------------------------------------------------------------
void PointCloudDefault::cleanup()
{
  // clear() keeps the capacity: shrink_to_fit() releases the memory (called after a failed allocation)
  release(Attr::All);
  coords.clear();
  coords.shrink_to_fit();
  n_points = 0;
  kdtree.reset();
}

// ------------------------------------------------------------
// Filter (in place)
// ------------------------------------------------------------
void PointCloudDefault::filter(const std::vector<bool>& keep)
{
  if (n_points != coords.size())
    throw std::logic_error("filter() cannot be called inside a ScopedPartition");
  if (keep.size() != n_points)
    throw std::invalid_argument("filter(): the mask must have one value per point");

  // Moves the kept values to the front (j <= i: nothing is overwritten before being read), then releases the end.
  auto compact = [&keep](auto& v)
  {
    if (v.empty()) return;
    size_t j = 0;
    for (size_t i = 0; i < v.size(); ++i)
    {
      if (keep[i]) v[j++] = v[i];
    }
    v.resize(j);
    v.shrink_to_fit();
  };

  compact(coords);
  visit_attributes([this, &compact](Attr flag, auto member, auto)
  {
    if (has(flag)) compact(this->*member);
  });

  n_points = coords.size();
  kdtree.reset();
}

// ------------------------------------------------------------
// Transforms
// ------------------------------------------------------------
void PointCloudDefault::translate(double x, double y, double z)
{
  for (size_t i = 0; i < n_points; ++i)
  {
    if (x != 0) coords[i].x -= static_cast<float>(x);
    if (y != 0) coords[i].y -= static_cast<float>(y);
    if (z != 0) coords[i].z -= static_cast<float>(z);
  }

  kdtree.reset();
}

void PointCloudDefault::scale(double x, double y, double z)
{
  for (size_t i = 0; i < n_points; ++i)
  {
    if (x != 1.0) coords[i].x *= static_cast<float>(x);
    if (y != 1.0) coords[i].y *= static_cast<float>(y);
    if (z != 1.0) coords[i].z *= static_cast<float>(z);
  }

  kdtree.reset();
}

// ------------------------------------------------------------
// Subset (by indices)
// ------------------------------------------------------------
PointCloudDefault PointCloudDefault::subset(const std::vector<unsigned int>& indices, Attr keep) const
{
  size_t new_count = indices.size();

  PointCloudDefault result;
  result.n_points = new_count;
  result.present = present & keep;

  auto gather = [&indices, new_count](auto& out, const auto& in)
  {
    out.resize(new_count);
    for (size_t j = 0; j < new_count; ++j)
      out[j] = in[indices[j]];
  };

  gather(result.coords, coords);
  visit_attributes([&](Attr flag, auto member, auto)
  {
    if (result.has(flag)) gather(result.*member, this->*member);
  });

  return result;
}

// ------------------------------------------------------------
// Merge operator (creates a new point cloud)
// ------------------------------------------------------------
PointCloudDefault PointCloudDefault::operator+(const PointCloudDefault& other) const
{
  PointCloudDefault result = *this;
  result += other;
  return result;
}

// ------------------------------------------------------------
// In-place merge operator
// ------------------------------------------------------------
PointCloudDefault& PointCloudDefault::operator+=(const PointCloudDefault& other)
{
  if (present != other.present)
    throw std::runtime_error("Cannot merge point clouds with different attributes: " + describe(present ^ other.present));

  return merge(other, present);
}

PointCloudDefault& PointCloudDefault::merge(const PointCloudDefault& other, Attr keep)
{
  if (n_points != coords.size() || other.n_points != other.coords.size())
    throw std::logic_error("merge() cannot be called inside a ScopedPartition");

  keep &= Attr::All;
  if (!has(keep))
    throw std::runtime_error("Cannot merge point clouds: missing attributes " + describe(keep & ~present));
  if (!other.has(keep))
    throw std::runtime_error("Cannot merge point clouds: the other point cloud is missing attributes " + describe(keep & ~other.present));

  // a += a: inserting a vector into itself is undefined
  if (&other == this)
  {
    const PointCloudDefault copy = other;
    return merge(copy, keep);
  }

  release(present & ~keep);

  coords.insert(coords.end(), other.coords.begin(), other.coords.end());
  visit_attributes([&](Attr flag, auto member, auto)
  {
    if (!has(flag)) return;
    auto& v = this->*member;
    const auto& o = other.*member;
    v.insert(v.end(), o.begin(), o.end());
  });

  n_points += other.n_points;
  kdtree.reset();

  return *this;
}

// ------------------------------------------------------------
// Swap (ScopedPartition, spatial_sort)
// ------------------------------------------------------------
void PointCloudDefault::swap_points(size_t i, size_t j)
{
  std::swap(coords[i], coords[j]);
  visit_attributes([&](Attr flag, auto member, auto)
  {
    if (has(flag)) std::swap((this->*member)[i], (this->*member)[j]);
  });
}

#endif

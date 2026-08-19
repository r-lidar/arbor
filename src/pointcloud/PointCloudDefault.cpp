#ifndef USING_R

#include "PointCloudDefault.h"

PointCloudDefault::PointCloudDefault()
{
  init();
}

PointCloudDefault::PointCloudDefault(size_t n, bool init_attributes)
{
  init();
  n_points = n;
  true_n_points = n;
  safe_alloc(n, init_attributes);
}


// ------------------------------------------------------------
// Cleanup
// ------------------------------------------------------------
void PointCloudDefault::cleanup()
{
  coords.clear();
  rgb.clear();
  treeid.clear();
  foliage.clear();
  passage.clear();
  hag.clear();
  pwood.clear();
  classif.clear();
  n_points = 0;
  true_n_points = 0;
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
PointCloudDefault PointCloudDefault::subset(const std::vector<unsigned int>& indices, bool xyz_only) const
{
  size_t new_count = indices.size();

  PointCloudDefault result;
  result.n_points = new_count;
  result.true_n_points = new_count;
  result.coords.resize(new_count);

  for (size_t j = 0; j < new_count; ++j)
  {
    unsigned int i = indices[j];
    result.coords[j] = coords[i];
  }

  if (!xyz_only)
  {
    if (!treeid.empty())
    {
      result.treeid.resize(new_count);
      for (size_t j = 0; j < new_count; ++j)
        result.treeid[j] = treeid[indices[j]];
    }

    if (!foliage.empty())
    {
      result.foliage.resize(new_count);
      for (size_t j = 0; j < new_count; ++j)
        result.foliage[j] = foliage[indices[j]];
    }

    if (!passage.empty())
    {
      result.passage.resize(new_count);
      for (size_t j = 0; j < new_count; ++j)
        result.passage[j] = passage[indices[j]];
    }

    if (!hag.empty())
    {
      result.hag.resize(new_count);
      for (size_t j = 0; j < new_count; ++j)
        result.hag[j] = hag[indices[j]];
    }

    if (!pwood.empty())
    {
      result.pwood.resize(new_count);
      for (size_t j = 0; j < new_count; ++j)
        result.pwood[j] = pwood[indices[j]];
    }

    if (!rgb.empty())
    {
      result.rgb.resize(new_count);
      for (size_t j = 0; j < new_count; ++j)
        result.rgb[j] = rgb[indices[j]];
    }

    if (!classif.empty())
    {
      result.classif.resize(new_count);
      for (size_t j = 0; j < new_count; ++j)
        result.classif[j] = classif[indices[j]];
    }

    if (!userdata.empty())
    {
        result.userdata.resize(new_count);
        for (size_t j = 0; j < new_count; ++j)
            result.userdata[j] = userdata[indices[j]];
    }
  }

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
  if (other.n_points == 0)
    return *this;

  size_t old_size = this->n_points;
  size_t new_size = old_size + other.n_points;

  // Merge coordinates
  coords.insert(coords.end(), other.coords.begin(), other.coords.end());

  // Merge rgb if both have it
  if (!rgb.empty() && !other.rgb.empty())
  {
    rgb.insert(rgb.end(), other.rgb.begin(), other.rgb.end());
  }
  else if (!rgb.empty())
  {
    rgb.clear();
  }

  // Merge treeid if both have it
  if (has_treeid() && other.has_treeid())
  {
    treeid.insert(treeid.end(), other.treeid.begin(), other.treeid.end());
  }
  else if (has_treeid())
  {
    treeid.clear();
  }

  // Merge foliage if both have it
  if (has_foliage() && other.has_foliage())
  {
    foliage.insert(foliage.end(), other.foliage.begin(), other.foliage.end());
  }
  else if (has_foliage())
  {
    foliage.clear();
  }

  // Merge passage if both have it
  if (has_passage() && other.has_passage())
  {
    passage.insert(passage.end(), other.passage.begin(), other.passage.end());
  }
  else if (has_passage())
  {
    passage.clear();
  }

  // Merge hag if both have it
  if (has_hag() && other.has_hag())
  {
    hag.insert(hag.end(), other.hag.begin(), other.hag.end());
  }
  else if (has_hag())
  {
    hag.clear();
  }

  // Merge pwood if both have it
  if (has_pwood() && other.has_pwood())
  {
    pwood.insert(pwood.end(), other.pwood.begin(), other.pwood.end());
  }
  else if (has_pwood())
  {
    pwood.clear();
  }

  // Merge classification if both have it
  if (has_class() && other.has_class())
  {
    classif.insert(classif.end(), other.classif.begin(), other.classif.end());
  }
  else if (has_class())
  {
    classif.clear();
  }

  // Merge classification if both have it
  if (has_userdata() && other.has_userdata())
  {
      userdata.insert(userdata.end(), other.userdata.begin(), other.userdata.end());
  }
  else if (has_userdata())
  {
      userdata.clear();
  }

  this->n_points = new_size;
  this->true_n_points = new_size;

  kdtree.reset();

  return *this;
}

// ------------------------------------------------------------
// Init
// ------------------------------------------------------------
void PointCloudDefault::init()
{
  n_points = 0;
  true_n_points = 0;
  coords.clear();
  rgb.clear();
  treeid.clear();
  foliage.clear();
  passage.clear();
  hag.clear();
  pwood.clear();
  classif.clear();
  userdata.clear();
}

// ------------------------------------------------------------
// Safe Alloc
// ------------------------------------------------------------
void PointCloudDefault::safe_alloc(size_t n, bool alloc_attrs)
{
  try
  {
    coords.resize(n, {0.0f, 0.0f, 0.0f});

    if (alloc_attrs)
    {
      treeid.resize(n, -1);
      foliage.resize(n, -1);
      classif.resize(n, 0);
      passage.resize(n, 0);
      hag.resize(n, 0.0f);
      pwood.resize(n, 0.0f);
      rgb.resize(n, {128, 128, 128});
      userdata.resize(n, 0);
    }
  }
  catch (...)
  {
    cleanup();
    throw;
  }
}

#endif

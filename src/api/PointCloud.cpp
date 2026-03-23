#include "PointCloud.h"

#include <algorithm> // for std::copy
#include <cstring>   // for std::memset
#include <algorithm>

#ifdef USING_R

PointCloudDataFrame::PointCloudDataFrame()
{
  init();
}

PointCloudDataFrame::PointCloudDataFrame(size_t n, bool init_attributes)
{
  init();
  n_points = n;
  true_n_points = n;
  safe_alloc(n, init_attributes);
}

PointCloudDataFrame::PointCloudDataFrame(const Rcpp::DataFrame& df)
{
  init();
  n_points = df.rows();
  true_n_points = n_points;
  owns_memory = false;

  std::vector<std::string> coord_names = {"X", "Y", "Z"};
  std::string treeid_name  = "treeID";
  std::string pwood_name   = "pwood";
  std::string foliage_name = "foliage";
  std::string hag_name     = "hag";
  std::string passage_name = "passage";
  std::string classif_name = "Classification";

  // --- Mandatory coordinates ---
  for (size_t i = 0; i < 3; ++i)
  {
    if (!df.containsElementNamed(coord_names[i].c_str()))
      throw std::runtime_error("Missing mandatory coordinate column: " + coord_names[i]);

    Rcpp::NumericVector col = df[coord_names[i]];
    coords[i] = col.begin();
  }

  // --- Optional attributes ---
  if (df.containsElementNamed(hag_name.c_str()))
  {
    Rcpp::NumericVector col = df[hag_name];
    hag = col.begin();
  }

  if (df.containsElementNamed(treeid_name.c_str()))
  {
    Rcpp::IntegerVector col = df[treeid_name];
    treeid = col.begin();
  }

  if (df.containsElementNamed(pwood_name.c_str()))
  {
    Rcpp::NumericVector col = df[pwood_name];
    pwood = col.begin();
  }

  if (df.containsElementNamed(foliage_name.c_str()))
  {
    Rcpp::IntegerVector col = df[foliage_name];
    foliage = col.begin();
  }

  if (df.containsElementNamed(passage_name.c_str()))
  {
    Rcpp::IntegerVector col = df[passage_name];
    passage = col.begin();
  }

  if (df.containsElementNamed(classif_name.c_str()))
  {
    Rcpp::IntegerVector col = df[classif_name];
    classif = col.begin();
  }
}

// ------------------------------------------------------------
// Destructor (Rule of Five: 1/5)
// ------------------------------------------------------------
PointCloudDataFrame::~PointCloudDataFrame()
{
  cleanup();
}

// ------------------------------------------------------------
// Copy Constructor (Rule of Five: 2/5)
// ------------------------------------------------------------
PointCloudDataFrame::PointCloudDataFrame(const PointCloudDataFrame& other)
{
  init();
  n_points = other.n_points;
  owns_memory = true; // A copy ALWAYS owns its new memory

  try
  {
    for (int i = 0; i < 3; ++i)
    {
      if (other.coords[i]) {
        coords[i] = new double[n_points];
        std::copy(other.coords[i], other.coords[i] + n_points, coords[i]);
      }
    }

    if (other.treeid) {
      treeid = new int[n_points];
      std::copy(other.treeid, other.treeid + n_points, treeid);
    }
    if (other.foliage) {
      foliage = new int[n_points];
      std::copy(other.foliage, other.foliage + n_points, foliage);
    }
    if (other.pwood) {
      pwood = new double[n_points];
      std::copy(other.pwood, other.pwood + n_points, pwood);
    }
    if (other.hag) {
      hag = new double[n_points];
      std::copy(other.hag, other.hag + n_points, hag);
    }
    if (other.passage) {
      passage = new int[n_points];
      std::copy(other.passage, other.passage + n_points, passage);
    }
    if (other.classif) {
      classif = new int[n_points];
      std::copy(other.classif, other.classif + n_points, classif);
    }
  }
  catch (...)
  {
    cleanup();
    throw;
  }
}

// ------------------------------------------------------------
// Move Constructor (Rule of Five: 3/5)
// ------------------------------------------------------------
PointCloudDataFrame::PointCloudDataFrame(PointCloudDataFrame&& other) noexcept
{
  init();
  swap(*this, other);
}

// ------------------------------------------------------------
// Copy Assignment Operator (Rule of Five: 4/5)
// ------------------------------------------------------------
// Using Copy-and-Swap.
// Note: The argument is passed by VALUE.
// If passed an l-value, Copy Constructor is called.
// If passed an r-value (temporary), Move Constructor is called.
PointCloudDataFrame& PointCloudDataFrame::operator=(PointCloudDataFrame other) noexcept
{
  swap(*this, other);
  return *this;
}

// ------------------------------------------------------------
// Move Assignment Operator (Rule of Five: 5/5) - OPTIONAL here
// ------------------------------------------------------------
// Because we use the "pass-by-value" idiom in operator= above,
// a specific Move Assignment (operator=(T&&)) is technically redundant
// but can be added for explicit clarity.
// The code above covers both cases.

// ------------------------------------------------------------
// Cleanup
// ------------------------------------------------------------
void PointCloudDataFrame::cleanup()
{
  if (!owns_memory)
    return;

  delete[] coords[0];
  delete[] coords[1];
  delete[] coords[2];

  delete[] treeid;
  delete[] foliage;
  delete[] pwood;
  delete[] hag;
  delete[] passage;
  delete[] classif;

  n_points = 0;
  true_n_points = 0;

  coords[0] = coords[1] = coords[2] = nullptr;
  treeid = foliage = classif = nullptr;
  pwood = hag = nullptr;
}

// ------------------------------------------------------------
// Transforms
// ------------------------------------------------------------
void PointCloudDataFrame::translate(double x, double y, double z)
{
  for (size_t i = 0; i < n_points; ++i)
  {
    if (x != 0) const_cast<double&>(coords[0][i]) -= x;
    if (y != 0) const_cast<double&>(coords[1][i]) -= y;
    if (z != 0) const_cast<double&>(coords[2][i]) -= z;
  }
}

void PointCloudDataFrame::scale(double x, double y, double z)
{
  for (size_t i = 0; i < n_points; ++i)
  {
    if (x != 1.0) const_cast<double&>(coords[0][i]) *= x;
    if (y != 1.0) const_cast<double&>(coords[1][i]) *= y;
    if (z != 1.0) const_cast<double&>(coords[2][i]) *= z;
  }
}

// ------------------------------------------------------------
// Subset
// ------------------------------------------------------------
// 1. New version using indices
PointCloudDataFrame PointCloudDataFrame::subset(const std::vector<int>& indices, bool xyz_only) const
{
  size_t new_count = indices.size();

  PointCloudDataFrame result;
  result.n_points = new_count;
  result.owns_memory = true;

  // Allocate Coords
  result.coords[0] = new double[new_count];
  result.coords[1] = new double[new_count];
  result.coords[2] = new double[new_count];

  // Allocate Optional Attributes
  if (!xyz_only)
  {
    if (treeid)  result.treeid  = new int[new_count];
    if (foliage) result.foliage = new int[new_count];
    if (passage) result.passage = new int[new_count];
    if (pwood)   result.pwood   = new double[new_count];
    if (hag)     result.hag     = new double[new_count];
  }

  for (size_t j = 0; j < new_count; ++j)
  {
    int i = indices[j];

    // XYZ
    result.coords[0][j] = coords[0][i];
    result.coords[1][j] = coords[1][i];
    result.coords[2][j] = coords[2][i];

    // Attributes
    if (!xyz_only)
    {
      if (treeid)  result.treeid[j]  = treeid[i];
      if (foliage) result.foliage[j] = foliage[i];
      if (passage) result.passage[j] = passage[i];
      if (pwood)   result.pwood[j]   = pwood[i];
      if (hag)     result.hag[j]     = hag[i];
    }
  }

  return result;
}

PointCloudDataFrame PointCloudDataFrame::subset(const std::vector<bool>& keep, bool xyz_only) const
{
  if (keep.size() != n_points)
    throw std::runtime_error("subset mask size mismatch: expected " + std::to_string(n_points) + " but got " + std::to_string(keep.size()));

  // Convert bool mask to indices
  std::vector<int> indices;
  indices.reserve(n_points/10);

  for (size_t i = 0; i < n_points; ++i) {
    if (keep[i]) {
      indices.push_back(i);
    }
  }

  // Handle the "all points" case early to avoid unnecessary copy/allocation
  if (indices.size() == n_points) {
    return *this;
  }

  return subset(indices, xyz_only);
}
// ------------------------------------------------------------
// Merge
// ------------------------------------------------------------

// Merge operator - creates a new point cloud by copying and using +=
PointCloudDataFrame PointCloudDataFrame::operator+(const PointCloudDataFrame& other) const
{
  PointCloudDataFrame result = *this;
  result += other;
  return result;
}

// In-place merge operator
PointCloudDataFrame& PointCloudDataFrame::operator+=(const PointCloudDataFrame& other)
{
  if (!this->owns_memory) {
    throw std::runtime_error("This point cloud does not own its memory because it is owned by R. Cannot merge.");
  }

  if (other.n_points == 0) {
    return *this;
  }

  size_t old_size = this->n_points;
  size_t new_size = old_size + other.n_points;

  // Reallocate and merge coordinates
  for (size_t d = 0; d < 3; ++d)
  {
    double* new_coords = new double[new_size];
    std::memcpy(new_coords, this->coords[d], old_size * sizeof(double));
    std::memcpy(new_coords + old_size, other.coords[d], other.n_points * sizeof(double));

    if (this->owns_memory)
      delete[] this->coords[d];

    this->coords[d] = new_coords;
  }

  // Merge treeid if both have it
  if (this->has_treeid() && other.has_treeid())
  {
    int* new_treeid = new int[new_size];
    std::memcpy(new_treeid, this->treeid, old_size * sizeof(int));
    std::memcpy(new_treeid + old_size, other.treeid, other.n_points * sizeof(int));

    if (this->owns_memory)
      delete[] this->treeid;

    this->treeid = new_treeid;
  }
  else if (this->has_treeid())
  {
    // This has treeid but other doesn't - remove treeid
    if (this->owns_memory) {
      delete[] this->treeid;
    }
    this->treeid = nullptr;
  }

  // Merge foliage if both have it
  if (this->has_foliage() && other.has_foliage())
  {
    int* new_foliage = new int[new_size];
    std::memcpy(new_foliage, this->foliage, old_size * sizeof(int));
    std::memcpy(new_foliage + old_size, other.foliage, other.n_points * sizeof(int));

    if (this->owns_memory)
      delete[] this->foliage;

    this->foliage = new_foliage;
  }
  else if (this->has_foliage())
  {
    if (this->owns_memory)
      delete[] this->foliage;

    this->foliage = nullptr;
  }

  // Merge passage if both have it
  if (this->has_passage() && other.has_passage())
  {
    int* new_passage = new int[new_size];
    std::memcpy(new_passage, this->passage, old_size * sizeof(int));
    std::memcpy(new_passage + old_size, other.passage, other.n_points * sizeof(int));

    if (this->owns_memory)
      delete[] this->passage;

    this->passage = new_passage;
  }
  else if (this->has_passage())
  {
    if (this->owns_memory)
      delete[] this->passage;

    this->passage = nullptr;
  }

  // Merge hag if both have it
  if (this->has_hag() && other.has_hag())
  {
    double* new_hag = new double[new_size];
    std::memcpy(new_hag, this->hag, old_size * sizeof(double));
    std::memcpy(new_hag + old_size, other.hag, other.n_points * sizeof(double));

    if (this->owns_memory)
      delete[] this->hag;

    this->hag = new_hag;
  }
  else if (this->has_hag())
  {
    if (this->owns_memory)
      delete[] this->hag;

    this->hag = nullptr;
  }

  // Merge pwood if both have it
  if (this->has_pwood() && other.has_pwood())
  {
    double* new_pwood = new double[new_size];
    std::memcpy(new_pwood, this->pwood, old_size * sizeof(double));
    std::memcpy(new_pwood + old_size, other.pwood, other.n_points * sizeof(double));

    if (this->owns_memory)
      delete[] this->pwood;

    this->pwood = new_pwood;
  }
  else if (this->has_pwood())
  {
    if (this->owns_memory)
      delete[] this->pwood;

    this->pwood = nullptr;
  }

  this->n_points = new_size;
  this->owns_memory = true;

  return *this;
}

void PointCloudDataFrame::swap(PointCloudDataFrame& first, PointCloudDataFrame& second) noexcept
{
  std::swap(first.n_points, second.n_points);
  std::swap(first.owns_memory, second.owns_memory);

  std::swap(first.coords[0], second.coords[0]);
  std::swap(first.coords[1], second.coords[1]);
  std::swap(first.coords[2], second.coords[2]);

  std::swap(first.treeid, second.treeid);
  std::swap(first.foliage, second.foliage);
  std::swap(first.passage, second.passage);
  std::swap(first.hag, second.hag);
  std::swap(first.pwood, second.pwood);
  std::swap(first.classif, second.classif);
}

void PointCloudDataFrame::init()
{
  n_points = 0;
  owns_memory = true;
  coords[0] = nullptr;
  coords[1] = nullptr;
  coords[2] = nullptr;
  treeid  = nullptr;
  foliage = nullptr;
  passage = nullptr;
  hag     = nullptr;
  pwood   = nullptr;
  classif = nullptr;
}

void PointCloudDataFrame::safe_alloc(size_t n, bool alloc_attrs)
{
  try
  {
    coords[0] = new double[n]();
    coords[1] = new double[n]();
    coords[2] = new double[n]();

    if (alloc_attrs)
    {
      treeid  = new int[n]();
      foliage = new int[n]();
      passage = new int[n]();
      classif = new int[n]();
      hag     = new double[n]();
      pwood   = new double[n]();
    }
  }
  catch (...)
  {
    // If ANY allocation fails, we must manually free what was already
    // allocated because the destructor won't run for a failed constructor.
    cleanup();
    throw;
  }
}

#else

PointCloudDefault::PointCloudDefault()
{
  init();
}

PointCloudDefault::PointCloudDefault(size_t n, bool init_attributes)
{
  init();
  n_points = n;
  safe_alloc(n, init_attributes);
}

PointCloudDefault::~PointCloudDefault()
{
  cleanup();
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
}

void PointCloudDefault::scale(double x, double y, double z)
{
  for (size_t i = 0; i < n_points; ++i)
  {
    if (x != 1.0) coords[i].x *= static_cast<float>(x);
    if (y != 1.0) coords[i].y *= static_cast<float>(y);
    if (z != 1.0) coords[i].z *= static_cast<float>(z);
  }
}

// ------------------------------------------------------------
// Subset (by indices)
// ------------------------------------------------------------
PointCloudDefault PointCloudDefault::subset(const std::vector<int>& indices, bool xyz_only) const
{
  size_t new_count = indices.size();

  PointCloudDefault result;
  result.n_points = new_count;
  result.coords.resize(new_count);

  for (size_t j = 0; j < new_count; ++j)
  {
    int i = indices[j];
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
  }

  return result;
}

// ------------------------------------------------------------
// Subset (by bool mask)
// ------------------------------------------------------------
PointCloudDefault PointCloudDefault::subset(const std::vector<bool>& keep, bool xyz_only) const
{
  if (keep.size() != n_points)
    throw std::runtime_error("subset mask size mismatch: expected " +
      std::to_string(n_points) + " but got " + std::to_string(keep.size()));

  // Convert bool mask to indices
  std::vector<int> indices;
  indices.reserve(n_points / 10);

  for (int i = 0; i < (int)n_points; ++i)
  {
    if (keep[i])
      indices.push_back(i);
  }

  // Handle the "all points" case early to avoid unnecessary copy/allocation
  if (indices.size() == n_points)
    return *this;

  return subset(indices, xyz_only);
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

  this->n_points = new_size;

  return *this;
}

// ------------------------------------------------------------
// Swap
// ------------------------------------------------------------
void PointCloudDefault::swap(PointCloudDefault& first, PointCloudDefault& second) noexcept
{
  std::swap(first.n_points, second.n_points);
  std::swap(first.coords, second.coords);
  std::swap(first.rgb, second.rgb);
  std::swap(first.treeid, second.treeid);
  std::swap(first.foliage, second.foliage);
  std::swap(first.passage, second.passage);
  std::swap(first.hag, second.hag);
  std::swap(first.pwood, second.pwood);
}

// ------------------------------------------------------------
// Init
// ------------------------------------------------------------
void PointCloudDefault::init()
{
  n_points = 0;
  coords.clear();
  rgb.clear();
  treeid.clear();
  foliage.clear();
  passage.clear();
  hag.clear();
  pwood.clear();
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
      treeid.resize(n, 0);
      foliage.resize(n, 0);
      passage.resize(n, 0);
      hag.resize(n, 0.0f);
      pwood.resize(n, 0.0f);
    }
  }
  catch (...)
  {
    cleanup();
    throw;
  }
}

#endif

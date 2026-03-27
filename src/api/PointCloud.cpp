#include "PointCloud.h"

#include <algorithm> // for std::copy
#include <cstring>   // for std::memset
#include <algorithm>
#include <random>

// Helper to convert HCL (Polar CIELAB) to RGB
// Based on standard conversion formulas (D65 illuminant)
static void hcl_to_rgb(float h, float c, float l, uint8_t* R, uint8_t* G, uint8_t* B)
{
  // 1. Convert HCL to CIELAB
  float h_rad = h * M_PI / 180.0f;
  float L = l;
  float a = std::cos(h_rad) * c;
  float b = std::sin(h_rad) * c;

  // 2. Convert CIELAB to XYZ
  auto f_inv = [](float t) {
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

  // The macro now takes a 'mandatory' boolean to handle X, Y, Z logic
  #define LOAD_ATTR_SAFE(expected_sexp, r_name, member, mandatory)         \
    if (df.containsElementNamed(r_name)) {                                 \
      SEXP col_sexp = df[r_name];                                          \
      if (TYPEOF(col_sexp) != expected_sexp) {                             \
        throw std::runtime_error(std::string("Column '") + r_name +        \
                                 "' has wrong type. Copy avoided.");       \
      }                                                                    \
      Rcpp::Vector<expected_sexp> col(col_sexp);                           \
      member = col.begin();                                                \
    } else if (mandatory) {                                                \
      throw std::runtime_error(std::string("Missing mandatory column: ")   \
                                 + r_name);                                \
    }

  // --- Mandatory ---
  LOAD_ATTR_SAFE(REALSXP, "X", coords[0], true)
  LOAD_ATTR_SAFE(REALSXP, "Y", coords[1], true)
  LOAD_ATTR_SAFE(REALSXP, "Z", coords[2], true)

  // --- Optional ---
  LOAD_ATTR_SAFE(REALSXP, "hag",            hag,     false)
  LOAD_ATTR_SAFE(INTSXP,  "treeID",         treeid,  false)
  LOAD_ATTR_SAFE(REALSXP, "pwood",          pwood,   false)
  LOAD_ATTR_SAFE(INTSXP,  "foliage",        foliage, false)
  LOAD_ATTR_SAFE(INTSXP,  "passage",        passage, false)
  LOAD_ATTR_SAFE(INTSXP,  "Classification", classif, false)
  LOAD_ATTR_SAFE(INTSXP,  "R",              red,     false)
  LOAD_ATTR_SAFE(INTSXP,  "G",              green,   false)
  LOAD_ATTR_SAFE(INTSXP,  "B",              blue,    false)

  #undef LOAD_ATTR_SAFE
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
    if (classif) result.classif = new int[new_count];
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
      if (classif) result.classif[j] = classif[i];
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
  if (!this->owns_memory)
  {
    throw std::runtime_error("This point cloud does not own its memory. Cannot merge.");
  }

  if (other.n_points == 0) return *this;

  size_t old_size = this->n_points;
  size_t other_size = other.n_points;

  // Merge coordinates (X, Y, Z)
  merge_attribute(this->coords[0], other.coords[0], old_size, other_size, true);
  merge_attribute(this->coords[1], other.coords[1], old_size, other_size, true);
  merge_attribute(this->coords[2], other.coords[2], old_size, other_size, true);

  // Merge optional attributes
  merge_attribute(this->treeid,  other.treeid,  old_size, other_size, other.has_treeid());
  merge_attribute(this->foliage, other.foliage, old_size, other_size, other.has_foliage());
  merge_attribute(this->passage, other.passage, old_size, other_size, other.has_passage());
  merge_attribute(this->hag,     other.hag,     old_size, other_size, other.has_hag());
  merge_attribute(this->pwood,   other.pwood,   old_size, other_size, other.has_pwood());
  merge_attribute(this->classif, other.classif, old_size, other_size, other.has_class());

  this->n_points += other_size;
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

void PointCloudDataFrame::colorize_trees(bool darken_foliage)
{
  if (!has_red() || !has_green() || !has_blue())
    throw std::runtime_error("RGB memory not allocated");

  if (!has_treeid())
    throw std::runtime_error("No treeID in this point cloud");

  if (!has_foliage())
    darken_foliage = false;

  if (true_size() == 0) return;

  struct RGB { uint8_t r,g,b; };

  std::unordered_map<int, RGB> color_cache;
  const float darken_factor = 0.7f;

  for (size_t i = 0; i < size(); ++i)
  {
    int id = get_treeid(i);
    if (id < 0) continue;

    auto [it, inserted] = color_cache.try_emplace(id, RGB{});
    if (inserted)
    {
      std::mt19937 gen(static_cast<uint32_t>(id));
      std::uniform_real_distribution<float> dist_h(0.0f, 360.0f);
      std::uniform_real_distribution<float> dist_c(42.0f, 98.0f);
      std::uniform_real_distribution<float> dist_l(40.0f, 90.0f);
      hcl_to_rgb(dist_h(gen), dist_c(gen), dist_l(gen), &it->second.r, &it->second.g, &it->second.b);
    }

    RGB color = it->second;

    if (darken_foliage && !is_wood(i))
    {
      color.r *= darken_factor;
      color.g *= darken_factor;
      color.b *= darken_factor;
    }
    set_red(i, static_cast<int>(color.r)*255);
    set_green(i, static_cast<int>(color.g)*255);
    set_blue(i, static_cast<int>(color.b)*255);
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
  result.true_n_points = new_count;
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

    if (!classif.empty())
    {
      result.classif.resize(new_count);
      for (size_t j = 0; j < new_count; ++j)
        result.classif[j] = classif[indices[j]];
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

  // Merge classification if both have it
  if (has_class() && other.has_class())
  {
    classif.insert(classif.end(), other.classif.begin(), other.classif.end());
  }
  else if (has_class())
  {
    classif.clear();
  }

  this->n_points = new_size;
  this->true_n_points = new_size;

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
    }
  }
  catch (...)
  {
    cleanup();
    throw;
  }
}

void PointCloudDefault::colorize_trees(bool darken_foliage)
{
    if (true_size() == 0) return;
    rgb.assign(true_size(), {170, 170, 170});

    std::unordered_map<int, RGB> color_cache;
    const float darken_factor = 0.7f;

    for (size_t i = 0; i < size(); ++i)
    {
        int id = get_treeid(i);
        if (id == -1) continue;

        auto [it, inserted] = color_cache.try_emplace(id, RGB{});
        if (inserted)
        {
            std::mt19937 gen(static_cast<uint32_t>(id));
            std::uniform_real_distribution<float> dist_h(0.0f, 360.0f);
            std::uniform_real_distribution<float> dist_c(42.0f, 98.0f);
            std::uniform_real_distribution<float> dist_l(40.0f, 90.0f);
            hcl_to_rgb(dist_h(gen), dist_c(gen), dist_l(gen), &it->second.r, &it->second.g, &it->second.b);
        }

        RGB color = it->second;

        if (darken_foliage && foliage[i] >= 1)
        {
            color.r *= darken_factor;
            color.g *= darken_factor;
            color.b *= darken_factor;
        }
        rgb[i] = color;
    }
}
#endif

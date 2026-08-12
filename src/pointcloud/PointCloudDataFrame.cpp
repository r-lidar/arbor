#ifdef USING_R

#include "PointCloudDataFrame.h"

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
  #define LOAD_ATTR_SAFE(expected_sexp, r_name, member, mandatory)       \
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
  true_n_points = other.true_n_points;
  owns_memory = true; // A copy ALWAYS owns its new memory
  kdtree = nullptr;

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
    if (other.red)   {
      red = new int[n_points];
      std::copy(other.red, other.red + n_points, red);
    }
    if (other.green) {
      green = new int[n_points];
      std::copy(other.green, other.green + n_points, green);
    }
    if (other.blue)  {
      blue  = new int[n_points];
      std::copy(other.blue, other.blue + n_points, blue);
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
: PointCloudDataFrame()
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
  delete[] red;
  delete[] green;
  delete[] blue;

  n_points = 0;
  true_n_points = 0;

  coords[0] = coords[1] = coords[2] = nullptr;
  treeid = foliage = classif = nullptr;
  pwood = hag = nullptr;
  red = green = blue = nullptr;

  kdtree.reset();
}

// ------------------------------------------------------------
// Transforms
// ------------------------------------------------------------
void PointCloudDataFrame::translate(double x, double y, double z)
{
  for (size_t i = 0; i < n_points; ++i)
  {
    if (x != 0) coords[0][i] -= x;
    if (y != 0) coords[1][i] -= y;
    if (z != 0) coords[2][i] -= z;
  }

  kdtree.reset();
}

void PointCloudDataFrame::scale(double x, double y, double z)
{
  for (size_t i = 0; i < n_points; ++i)
  {
    if (x != 1.0) coords[0][i] *= x;
    if (y != 1.0) coords[1][i] *= y;
    if (z != 1.0) coords[2][i] *= z;
  }

  kdtree.reset();
}

// ------------------------------------------------------------
// Subset
// ------------------------------------------------------------
PointCloudDataFrame PointCloudDataFrame::subset(const std::vector<unsigned int>& indices, bool xyz_only) const
{
  size_t new_count = indices.size();

  PointCloudDataFrame result;
  result.n_points = new_count;
  result.true_n_points = new_count;
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
    if (red)     result.red     = new int[new_count];
    if (green)   result.green   = new int[new_count];
    if (blue)    result.blue    = new int[new_count];
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
      if (red)     result.red[j]     = red[i];
      if (green)   result.green[j]   = green[i];
      if (blue)    result.blue [j]   = blue[i];
    }
  }

  return result;
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
  merge_attribute(this->red,     other.red,     old_size, other_size, other.has_red());
  merge_attribute(this->green,   other.green,   old_size, other_size, other.has_green());
  merge_attribute(this->blue,    other.blue,    old_size, other_size, other.has_blue());

  this->n_points += other_size;
  this->owns_memory = true;

  kdtree.reset();

  return *this;
}

void PointCloudDataFrame::swap(PointCloudDataFrame& first, PointCloudDataFrame& second) noexcept
{
  std::swap(first.n_points, second.n_points);
  std::swap(first.true_n_points, second.true_n_points);
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

  std::swap(first.red,   second.red);
  std::swap(first.green, second.green);
  std::swap(first.blue,  second.blue);
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

#endif

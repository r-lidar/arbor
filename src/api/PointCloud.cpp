#include "PointCloud.h"

#include <algorithm> // for std::copy
#include <cstring>   // for std::memset
#include <algorithm>

#ifdef USING_R

PointCloudDataFrame::PointCloudDataFrame(size_t n, bool init_attributes)
{
  if (n < 0) throw std::invalid_argument("PointCloudDataFrame: n must be >= 0");

  n_points = n;
  owns_memory = true;

  // Allocate coordinates
  for (int d = 0; d < 3; ++d)
    coords[d] = new double[n_points]();

    // Optional attributes
    if (init_attributes)
    {
      treeid  = new int[n_points]();
      foliage = new int[n_points]();
      passage = new int[n_points]();
      hag     = new double[n_points]();
      pwood   = new double[n_points]();
    }
}

PointCloudDataFrame::PointCloudDataFrame(const Rcpp::DataFrame& df)
{
  owns_memory = false;

  std::vector<std::string> coord_names = {"X", "Y", "Z"};
  std::string treeid_name  = "treeID";
  std::string pwood_name   = "pwood";
  std::string foliage_name = "foliage";
  std::string hag_name     = "hag";
  std::string passage_name = "passage";

  n_points = df.rows();

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
}

// ------------------------------------------------------------
// Copy constructor (deep copy)
// ------------------------------------------------------------
PointCloudDataFrame::PointCloudDataFrame(const PointCloudDataFrame& other)
{
  n_points = other.n_points;
  owns_memory = true;

  coords[0] = new double[n_points];
  coords[1] = new double[n_points];
  coords[2] = new double[n_points];

  std::copy(other.coords[0], other.coords[0] + n_points, coords[0]);
  std::copy(other.coords[1], other.coords[1] + n_points, coords[1]);
  std::copy(other.coords[2], other.coords[2] + n_points, coords[2]);

  if (other.treeid)
  {
    treeid = new int[n_points];
    std::copy(other.treeid, other.treeid + n_points, treeid);
  }

  if (other.foliage)
  {
    foliage = new int[n_points];
    std::copy(other.foliage, other.foliage + n_points, foliage);
  }

  if (other.pwood)
  {
    pwood = new double[n_points];
    std::copy(other.pwood, other.pwood + n_points, pwood);
  }

  if (other.hag)
  {
    hag = new double[n_points];
    std::copy(other.hag, other.hag + n_points, hag);
  }

  if (other.hag)
  {
    passage = new int[n_points];
    std::copy(other.passage, other.passage + n_points, passage);
  }
}

// ------------------------------------------------------------
// Move constructor
// ------------------------------------------------------------
PointCloudDataFrame::PointCloudDataFrame(PointCloudDataFrame&& other) noexcept
{
  n_points    = other.n_points;
  owns_memory = other.owns_memory;

  coords[0] = other.coords[0];
  coords[1] = other.coords[1];
  coords[2] = other.coords[2];
  treeid    = other.treeid;
  foliage   = other.foliage;
  pwood     = other.pwood;
  hag       = other.hag;
  passage   = other.passage;

  other.coords[0] = other.coords[1] = other.coords[2] = nullptr;
  other.treeid = other.foliage = nullptr;
  other.pwood = nullptr;
  other.hag = nullptr;
  other.passage = nullptr;
  other.owns_memory = false;
}

// ------------------------------------------------------------
// Copy assignment
// ------------------------------------------------------------
PointCloudDataFrame& PointCloudDataFrame::operator=(const PointCloudDataFrame& other)
{
  if (this != &other)
  {
    if (owns_memory) cleanup();

    *this = PointCloudDataFrame(other);
  }
  return *this;
}

// ------------------------------------------------------------
// Move assignment
// ------------------------------------------------------------
PointCloudDataFrame& PointCloudDataFrame::operator=(PointCloudDataFrame&& other) noexcept
{
  if (this != &other)
  {
    if (owns_memory) cleanup();

    n_points    = other.n_points;
    owns_memory = other.owns_memory;

    coords[0] = other.coords[0];
    coords[1] = other.coords[1];
    coords[2] = other.coords[2];
    treeid    = other.treeid;
    foliage   = other.foliage;
    pwood     = other.pwood;
    hag       = other.hag;
    passage   = other.passage;

    other.coords[0] = other.coords[1] = other.coords[2] = nullptr;
    other.treeid = other.foliage = nullptr;
    other.pwood = nullptr;
    other.hag = nullptr;
    other.passage = nullptr;
    other.owns_memory = false;
  }
  return *this;
}

// ------------------------------------------------------------
// Destructor
// ------------------------------------------------------------
PointCloudDataFrame::~PointCloudDataFrame()
{
  cleanup();
}

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

  coords[0] = coords[1] = coords[2] = nullptr;
  treeid = foliage = nullptr;
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
PointCloudDataFrame PointCloudDataFrame::subset(const std::vector<bool>& keep, bool xyz_only) const
{
  if (keep.size() != n_points)
    throw std::runtime_error("subset mask size mismatch: expected " + std::to_string(n_points) + " but got " + std::to_string(keep.size()));

  size_t new_count = std::count(keep.begin(), keep.end(), true);

  if (new_count == n_points)
    return *this;

  PointCloudDataFrame result;
  result.n_points = new_count;
  result.owns_memory = true;

  result.coords[0] = new double[new_count];
  result.coords[1] = new double[new_count];
  result.coords[2] = new double[new_count];

  size_t j = 0;
  for (size_t i = 0; i < n_points; ++i)
  {
    if (keep[i])
    {
      result.coords[0][j] = coords[0][i];
      result.coords[1][j] = coords[1][i];
      result.coords[2][j] = coords[2][i];
      ++j;
    }
  }

  if (xyz_only)  return result;

  if (treeid)  result.treeid  = new int[new_count];
  if (foliage) result.foliage = new int[new_count];
  if (passage) result.passage = new int[new_count];
  if (pwood)   result.pwood   = new double[new_count];
  if (hag)     result.hag     = new double[new_count];


  j = 0;
  for (size_t i = 0; i < n_points; ++i)
  {
    if (keep[i])
    {
      if (treeid)  result.treeid[j]  = treeid[i];
      if (foliage) result.foliage[j] = foliage[i];
      if (passage) result.passage[j] = passage[i];
      if (pwood)   result.pwood[j]   = pwood[i];
      if (hag)     result.hag[j]     = hag[i];
      ++j;
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
#endif

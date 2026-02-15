#include "Adaptor.h"

#ifdef USING_R
PointCloud::PointCloud(const Rcpp::DataFrame& df)
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
PointCloud::PointCloud(const PointCloud& other)
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
PointCloud::PointCloud(PointCloud&& other) noexcept
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
PointCloud& PointCloud::operator=(const PointCloud& other)
{
  if (this != &other)
  {
    if (owns_memory) cleanup();

    *this = PointCloud(other);
  }
  return *this;
}

// ------------------------------------------------------------
// Move assignment
// ------------------------------------------------------------
PointCloud& PointCloud::operator=(PointCloud&& other) noexcept
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
PointCloud::~PointCloud()
{
  cleanup();
}

// ------------------------------------------------------------
// Cleanup
// ------------------------------------------------------------
void PointCloud::cleanup()
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
void PointCloud::translate(double x, double y, double z)
{
  for (size_t i = 0; i < n_points; ++i)
  {
    if (x != 0) const_cast<double&>(coords[0][i]) -= x;
    if (y != 0) const_cast<double&>(coords[1][i]) -= y;
    if (z != 0) const_cast<double&>(coords[2][i]) -= z;
  }
}

void PointCloud::scale(double x, double y, double z)
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
PointCloud PointCloud::subset(const std::vector<bool>& keep, bool xyz_only) const
{
  if (keep.size() != n_points)
    throw std::runtime_error("subset mask size mismatch: expected " + std::to_string(n_points) + " but got " + std::to_string(keep.size()));

  size_t new_count = std::count(keep.begin(), keep.end(), true);

  if (new_count == n_points)
    return *this;

  PointCloud result;
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
#endif

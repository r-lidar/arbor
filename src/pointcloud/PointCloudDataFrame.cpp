#ifdef USING_R

#include "PointCloudDataFrame.h"

#include <type_traits>

// R type of a column from the type of its member: int* is INTSXP, double* is REALSXP
template <typename T> constexpr int sexp_type_of() { return std::is_same_v<T, int> ? INTSXP : REALSXP; }

PointCloudDataFrame::PointCloudDataFrame()
{
}

PointCloudDataFrame::PointCloudDataFrame(size_t n, Attr attrs)
{
  n_points = n;
  safe_alloc(n, attrs);
}

PointCloudDataFrame::PointCloudDataFrame(const Rcpp::DataFrame& df)
{
  n_points = df.rows();
  owns_memory = false;

  // swap_points() relies on every column being REALSXP, INTSXP, or LGLSXP.
  // Fail fast at construction rather than silently desyncing a column
  // (e.g. a STRSXP or list column) during a ScopedPartition (swap_points) later on.
  Rcpp::CharacterVector names = df.names();
  for (int c = 0; c < df.size(); ++c)
  {
    SEXP col = df[c];
    switch (TYPEOF(col))
    {
    case REALSXP:
    case INTSXP:
    case LGLSXP:
      break;
    default:
      throw std::runtime_error(std::string("Column '") + Rcpp::as<std::string>(names[c]) + "' has unsupported type for swapping " + "(only numeric, integer, and logical columns are supported).");
    }
  }

  // This is for sorting R data.frames
  generic_columns.reserve(df.size());
  for (int c = 0; c < df.size(); ++c)
  {
    SEXP col = df[c];
    switch (TYPEOF(col))
    {
    case REALSXP: generic_columns.push_back({REAL(col),    REALSXP}); break;
    case INTSXP:  generic_columns.push_back({INTEGER(col), INTSXP});  break;
    case LGLSXP:  generic_columns.push_back({LOGICAL(col), LGLSXP});  break;
    }
  }

  // --- Mandatory ---
  const char* xyz[3] = {"X", "Y", "Z"};
  for (int d = 0; d < 3; ++d)
  {
    if (!df.containsElementNamed(xyz[d]))
      throw std::runtime_error(std::string("Missing mandatory column: ") + xyz[d]);
    SEXP col = df[xyz[d]];
    if (TYPEOF(col) != REALSXP)
      throw std::runtime_error(std::string("Column '") + xyz[d] + "' has wrong type. Copy avoided.");
    coords[d] = REAL(col);
  }

  // --- Optional ---
  // An attribute is present only if all its columns are (RGB is three columns)
  Attr missing = Attr::None;
  visit_attributes([&](Attr flag, auto member, auto, const char* r_name)
  {
    using T = std::remove_pointer_t<std::remove_reference_t<decltype(this->*member)>>;
    if (!df.containsElementNamed(r_name))
    {
      missing |= flag;
      return;
    }

    // Not used by the pipeline: a column of another type is ignored rather than rejected
    constexpr Attr optional = Attr::Intensity | Attr::GPStime | Attr::ReturnNumber | Attr::NumberOfReturns;

    SEXP col = df[r_name];
    if (TYPEOF(col) != sexp_type_of<T>() && any(flag & optional))
    {
      missing |= flag;
      return;
    }
    if (TYPEOF(col) != sexp_type_of<T>())
      throw std::runtime_error(std::string("Column '") + r_name + "' has wrong type. Copy avoided.");

    if constexpr (std::is_same_v<T, int>) this->*member = INTEGER(col);
    else this->*member = REAL(col);
    present |= flag;
  });
  present &= ~missing;
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
  n_points = other.n_points;
  owns_memory = true; // A copy ALWAYS owns its new memory

  try
  {
    auto copy = [this](auto*& out, const auto* in)
    {
      using T = std::remove_const_t<std::remove_pointer_t<decltype(in)>>;
      out = new T[n_points];
      std::copy(in, in + n_points, out);
    };

    for (int i = 0; i < 3; ++i) copy(coords[i], other.coords[i]);

    visit_attributes([&](Attr flag, auto member, auto, const char*)
    {
      if (other.has(flag)) copy(this->*member, other.*member);
    });
    present = other.present;
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
  // Every non-null pointer is owned when owns_memory: the columns of a failed allocation included
  for (int d = 0; d < 3; ++d)
  {
    if (owns_memory) delete[] coords[d];
    coords[d] = nullptr;
  }

  visit_attributes([this](Attr, auto member, auto, const char*)
  {
    if (owns_memory) delete[] (this->*member);
    this->*member = nullptr;
  });

  n_points = 0;
  present = Attr::None;
  owns_memory = true;
  generic_columns.clear();

  kdtree.reset();
}

// ------------------------------------------------------------
// Attribute management
// ------------------------------------------------------------
void PointCloudDataFrame::allocate(Attr)
{
  throw std::runtime_error("PointCloudDataFrame cannot allocate attributes after construction");
}

void PointCloudDataFrame::release(Attr)
{
  throw std::runtime_error("PointCloudDataFrame cannot release attributes after construction");
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
PointCloudDataFrame PointCloudDataFrame::subset(const std::vector<unsigned int>& indices, Attr keep) const
{
  size_t new_count = indices.size();

  PointCloudDataFrame result;
  result.n_points = new_count;
  result.owns_memory = true;

  auto gather = [&indices, new_count](auto*& out, const auto* in)
  {
    using T = std::remove_const_t<std::remove_pointer_t<decltype(in)>>;
    out = new T[new_count];
    for (size_t j = 0; j < new_count; ++j)
      out[j] = in[indices[j]];
  };

  // result owns what is already allocated: its destructor frees it if an allocation throws
  for (int d = 0; d < 3; ++d) gather(result.coords[d], coords[d]);

  visit_attributes([&](Attr flag, auto member, auto, const char*)
  {
    if (has(flag) && any(keep & flag)) gather(result.*member, this->*member);
  });
  result.present = present & keep;

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
  if (present != other.present)
    throw std::runtime_error("Cannot merge point clouds with different attributes: " + describe(present ^ other.present));

  return merge(other, present);
}

PointCloudDataFrame& PointCloudDataFrame::merge(const PointCloudDataFrame& other, Attr keep)
{
  if (!owns_memory)
    throw std::runtime_error("This point cloud does not own its memory. Cannot merge.");

  keep &= Attr::All;
  if (!has(keep))
    throw std::runtime_error("Cannot merge point clouds: missing attributes " + describe(keep & ~present));
  if (!other.has(keep))
    throw std::runtime_error("Cannot merge point clouds: the other point cloud is missing attributes " + describe(keep & ~other.present));

  const size_t old_n = n_points;
  const size_t other_n = other.n_points;

  // Reallocates and concatenates. Reads other's column before freeing this one: other may be *this.
  auto concat = [old_n, other_n](auto*& current, const auto* in)
  {
    using T = std::remove_const_t<std::remove_pointer_t<decltype(in)>>;
    T* next = new T[old_n + other_n];
    std::copy(current, current + old_n, next);
    std::copy(in, in + other_n, next + old_n);
    delete[] current;
    current = next;
  };

  // Allocate first, free the dropped columns last: other may be *this
  for (int d = 0; d < 3; ++d) concat(coords[d], other.coords[d]);
  visit_attributes([&](Attr flag, auto member, auto, const char*)
  {
    if (any(keep & flag)) concat(this->*member, other.*member);
  });
  visit_attributes([&](Attr flag, auto member, auto, const char*)
  {
    if (any(keep & flag)) return;
    delete[] (this->*member);
    this->*member = nullptr;
  });

  present = keep;
  n_points = old_n + other_n;
  kdtree.reset();

  return *this;
}

void PointCloudDataFrame::swap_points(size_t i, size_t j)
{
  if (!owns_memory)
  {
    // Backed directly by an R data.frame. Swap every column generically
    // so nothing gets left behind / desynced, regardless of whether this
    // class explicitly models that column.
    for (auto& gc : generic_columns)
    {
      switch (gc.type)
      {
      case REALSXP: std::swap(static_cast<double*>(gc.ptr)[i], static_cast<double*>(gc.ptr)[j]); break;
      case INTSXP:
      case LGLSXP:  std::swap(static_cast<int*>(gc.ptr)[i],    static_cast<int*>(gc.ptr)[j]);    break;
      default: break;
      }
    }
  }
  else
  {
    // Self-owned memory: no source_df, only the known attributes exist.
    for (int d = 0; d < 3; ++d) std::swap(coords[d][i], coords[d][j]);
    visit_attributes([&](Attr flag, auto member, auto, const char*)
    {
      if (has(flag)) std::swap((this->*member)[i], (this->*member)[j]);
    });
  }
}

void PointCloudDataFrame::swap(PointCloudDataFrame& first, PointCloudDataFrame& second) noexcept
{
  std::swap(first.n_points, second.n_points);
  std::swap(first.owns_memory, second.owns_memory);
  std::swap(first.present, second.present);

  std::swap(first.coords[0], second.coords[0]);
  std::swap(first.coords[1], second.coords[1]);
  std::swap(first.coords[2], second.coords[2]);

  visit_attributes([&](Attr, auto member, auto, const char*)
  {
    std::swap(first.*member, second.*member);
  });

  std::swap(first.generic_columns, second.generic_columns);
}

void PointCloudDataFrame::safe_alloc(size_t n, Attr attrs)
{
  try
  {
    for (int d = 0; d < 3; ++d) coords[d] = new double[n]();

    visit_attributes([&](Attr flag, auto member, auto def, const char*)
    {
      if (!any(attrs & flag)) return;
      using T = std::remove_pointer_t<std::remove_reference_t<decltype(this->*member)>>;
      this->*member = new T[n];
      std::fill(this->*member, this->*member + n, static_cast<T>(def));
      present |= flag;
    });
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

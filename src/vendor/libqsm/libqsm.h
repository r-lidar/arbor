/**
 * @file libqsm.h
 * Project: Arbor
 *
 * MIT License
 *
 * Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*
 *  Read usage
 *
 *   libqsm::QSMreader r("tree.qsm");
 *   for (const auto& node : r.nodes()) { ... node ... }
 *   for (const auto& edge : r.edges()) { ... edge ... }
 *
 *  Write usage
 *
 *   libqsm::QSMwriter w("out.qsm");
 *   w.set_software("MyApp 1.0");
 *   w.set_origin(448231.0, 5412087.0, 312.5);
 *   w.add_nodes(qn.begin(), qn.end());
 *   w.add_edges(qe.begin(), qe.end());
 *   w.write();
 *
 */

#ifndef LIBQSM_H
#define LIBQSM_H

#include <cstdint>
#include <string>
#include <vector>
#include <stdexcept>
#include <iosfwd>    // forward-declares std::ifstream / std::ofstream
#include <ostream>

#define LIBQSM_VERSION_MAJOR 1
#define LIBQSM_VERSION_MINOR 0

#define LIBQSM_STRINGIFY2(x) #x
#define LIBQSM_STRINGIFY(x) LIBQSM_STRINGIFY2(x)

#define LIBQSM_VERSION_STRING              \
LIBQSM_STRINGIFY(LIBQSM_VERSION_MAJOR) "." \
LIBQSM_STRINGIFY(LIBQSM_VERSION_MINOR)

namespace libqsm {

// ---------------------------------------------------------------------------
// Sentinel values
// ---------------------------------------------------------------------------

static constexpr float   UNSET_FLOAT = -1.0f;

// ---------------------------------------------------------------------------
// In-memory record structs
//
// These structs are the canonical in-memory representation and carry the
// UNION of all fields across every spec version.
//
// Fields introduced in versions later than 1.0 carry a sentinel default so
// that code reading a v1.0 file always sees a well-defined, detectable value.
//
// When adding a new field, annotate it with the minimum version tag:
//   float my_new_field = UNSET_FLOAT;   // [v1.1+]
// ---------------------------------------------------------------------------

struct QSMnode
{
  double  x  = 0.0;
  double  y  = 0.0;
  double  z  = 0.0;
};

struct QSMedge
{
  // v1.0
  uint32_t source = 0;
  uint32_t target = 0;
  float    radius = UNSET_FLOAT;
  uint8_t quality = 0;

  // v1.1
  float subtree_length   = UNSET_FLOAT;
  float distance_to_root = UNSET_FLOAT;

  // v1.2
  int32_t axis_id      = 0;
  uint8_t branch_order = 0;
};

// ---------------------------------------------------------------------------
// Header
// ---------------------------------------------------------------------------

struct QSMheader
{
  uint8_t  version_major = LIBQSM_VERSION_MAJOR;
  uint8_t  version_minor = LIBQSM_VERSION_MINOR;
  std::string software;
  std::string created;
  std::string crs;
  uint64_t    n_nodes   = 0;
  uint64_t    n_edges   = 0;
  uint32_t    node_size = 16;
  uint32_t    edge_size = 30;
  double      x_offset  = 0.0;
  double      y_offset  = 0.0;
  double      z_offset  = 0.0;
  std::vector<std::string> messages;
};


inline std::ostream& operator<<(std::ostream& os, const QSMnode& n)
{
  os << "{x=" << n.x
     << ", y=" << n.y
     << ", z=" << n.z
     << "}";
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const QSMedge& e)
{
  os << "{ source=" << e.source
     << ", target=" << e.target
     << ", radius=" << e.radius
     << ", subtree_length=" << e.subtree_length
     << ", distance_to_root=" << e.distance_to_root
     << ", axis_id=" << e.axis_id
     << ", branch_order=" << static_cast<int>(e.branch_order)
     << ", quality=" << static_cast<int>(e.quality)
     << "}";
  return os;
}
// ---------------------------------------------------------------------------
// Version dispatch table
//
// Each spec version owns four functions that speak *only* that version's
// binary layout.
//
//   read_*  reads exactly min_*_size mandatory bytes from the stream and
//           returns a fully-populated in-memory record; absent fields (those
//           not yet defined in that version) are left at their sentinel.
//
//   write_* writes exactly min_*_size bytes to the stream.
//
// The caller (read_binary / write_nodes / write_edges) is responsible for
// skipping any extra bytes beyond min_*_size that the file may carry (forward
// compatibility with a newer writer's extension fields).
// ---------------------------------------------------------------------------

struct VersionDispatch
{
  using NodeReader = QSMnode (*)(std::ifstream&,  const QSMheader&);
  using EdgeReader = QSMedge (*)(std::ifstream&,  const QSMheader&);
  using NodeWriter = void    (*)(std::ofstream&,  const QSMnode&,  const QSMheader&);
  using EdgeWriter = void    (*)(std::ofstream&,  const QSMedge&,  const QSMheader&);

  NodeReader read_node  = nullptr;
  EdgeReader read_edge  = nullptr;
  NodeWriter write_node = nullptr;
  EdgeWriter write_edge = nullptr;
};

struct QSMversionSpec
{
  uint8_t  major;
  uint8_t  minor;
  uint32_t min_node_size;    // mandatory bytes per node record for this version
  uint32_t min_edge_size;    // mandatory bytes per edge record for this version
  VersionDispatch dispatch;  // reader/writer functors for this version
};

// ---------------------------------------------------------------------------
// Registry and resolver - defined in libqsm.cpp
//
// version_registry() returns the ordered table of known versions (lowest
// first).  Adding a new version means appending one row here; nothing else
// changes in the library core.
//
// resolve_version() selects the best-matching entry:
//   - exact match               -> return that entry
//   - known major/unknown minor -> return the highest known minor, emit warning
//   - unknown major             -> throw (breaking change; cannot interpret safely)
// ---------------------------------------------------------------------------

const std::vector<QSMversionSpec>& version_registry();
const QSMversionSpec& resolve_version(uint8_t major, uint8_t minor, std::string* warning_out = nullptr);

// ---------------------------------------------------------------------------
// QSMreader
// ---------------------------------------------------------------------------

class QSMreader
{
public:
  explicit QSMreader(const std::string& filename);

  // Header metadata
  uint8_t            get_version_major() const noexcept { return header.version_major; }
  uint8_t            get_version_minor() const noexcept { return header.version_minor; }
  const std::string& get_software()      const noexcept { return header.software;  }
  const std::string& get_created()       const noexcept { return header.created;   }
  const std::string& get_crs()           const noexcept { return header.crs;       }
  double             get_x_offset()      const noexcept { return header.x_offset;  }
  double             get_y_offset()      const noexcept { return header.y_offset;  }
  double             get_z_offset()      const noexcept { return header.z_offset;  }
  uint32_t           get_node_size()     const noexcept { return header.node_size; }
  uint32_t           get_edge_size()     const noexcept { return header.edge_size; }
  int                get_message_count() const noexcept { return static_cast<int>(header.messages.size()); }
  const std::string& get_message(int i)  const          { return header.messages.at(static_cast<std::size_t>(i)); }

  // Record counts
  int node_count() const noexcept { return static_cast<int>(nodes_.size()); }
  int edge_count() const noexcept { return static_cast<int>(edges_.size()); }

  // Index access
  const QSMnode& node(int i) const { return nodes_.at(static_cast<std::size_t>(i)); }
  const QSMedge& edge(int i) const { return edges_.at(static_cast<std::size_t>(i)); }

  // Vector access (for range-for)
  const std::vector<QSMnode>& nodes() const noexcept { return nodes_; }
  const std::vector<QSMedge>& edges() const noexcept { return edges_; }

private:
  void parse();

  std::string    filename_;
  std::streampos data_start_;

  QSMheader            header;
  std::vector<QSMnode> nodes_;
  std::vector<QSMedge> edges_;
};

// ---------------------------------------------------------------------------
// QSMwriter
// ---------------------------------------------------------------------------

class QSMwriter
{
public:
  explicit QSMwriter(const std::string& filename);

  void set_version_major(uint8_t v) noexcept             { header.version_major = v; }
  void set_version_minor(uint8_t v) noexcept             { header.version_minor = v; }
  void set_software     (const std::string& s) noexcept  { header.software = s; }
  void set_crs          (const std::string& c) noexcept  { header.crs      = c; }
  void add_message      (const std::string& m) noexcept  { header.messages.push_back(m); }
  void set_origin(double x, double y, double z) noexcept { header.x_offset = x; header.y_offset = y; header.z_offset = z; }

  void add_node(const QSMnode& n) { nodes_.push_back(n); }
  void add_edge(const QSMedge& e) { edges_.push_back(e); }

  template <typename InputIt>
  void add_nodes(InputIt first, InputIt last) { nodes_.insert(nodes_.end(), first, last); }

  template <typename InputIt>
  void add_edges(InputIt first, InputIt last) { edges_.insert(edges_.end(), first, last); }

  void clear_records() noexcept { nodes_.clear(); edges_.clear(); }

  void write();

private:
  void write_header(std::ofstream& out, const QSMversionSpec& spec) const;
  void write_nodes (std::ofstream& out, const QSMversionSpec& spec) const;
  void write_edges (std::ofstream& out, const QSMversionSpec& spec) const;

  std::string filename_;

  QSMheader            header;
  std::vector<QSMnode> nodes_;
  std::vector<QSMedge> edges_;
};

} // namespace libqsm

#endif // LIBQSM_H

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
 *   w.set_format(2);           // 0 = minimal, 1 = + architecture, 2 = + distances
 *   w.add_node_extra_field(libqsm::ExtraFieldType::F32, "wood_density");
 *   w.add_edge_extra_field(libqsm::ExtraFieldType::F32, "taper_rate");
 *   w.add_nodes(qn.begin(), qn.end());
 *   w.add_edges(qe.begin(), qe.end());
 *   w.write();
 *   // Bounding box is computed automatically from nodes during write().
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
#include <unordered_map>

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

static constexpr float UNSET_FLOAT = -1.0f;

// ---------------------------------------------------------------------------
// In-memory record structs
//
// These structs are the canonical in-memory representation and carry the
// UNION of all fields across every format level.
//
// Fields introduced in format levels higher than 0 carry a sentinel default
// so that code reading a Format 0 file always sees a well-defined, detectable
// value.
//
// When adding a new field, annotate it with the minimum format tag:
//   float my_new_field = UNSET_FLOAT;   // [format 2+]
// ---------------------------------------------------------------------------

struct QSMnode
{
  double  x  = 0.0;
  double  y  = 0.0;
  double  z  = 0.0;
};

struct QSMedge
{
  // format 0  (offsets 0-12, 13 bytes total)
  uint32_t source  = 0;
  uint32_t target  = 0;
  float    radius  = UNSET_FLOAT;
  uint8_t  quality = 0;

  // format 1+  (offsets 13-17, +5 bytes)
  uint32_t axis_id      = 0;
  uint8_t  branch_order = 0;

  // format 2+  (offsets 18-25, +8 bytes)
  float subtree_length   = UNSET_FLOAT;
  float distance_to_root = UNSET_FLOAT;
};

// ---------------------------------------------------------------------------
// Extra-field type tokens
// ---------------------------------------------------------------------------

enum class ExtraFieldType : uint8_t
{
  I8, U8, I16, U16, I32, U32, F32, I64, U64, F64
};

// Returns the byte size of the type token.
uint32_t extra_field_type_size(ExtraFieldType t) noexcept;

// Returns the canonical string token for a type (e.g. ExtraFieldType::F32 -> "F32").
const char* extra_field_type_name(ExtraFieldType t) noexcept;

// Parses a value string of the form "<TYPE> <name>" (e.g. "F32 wood_density").
// Throws std::runtime_error on malformed input or unrecognised type token.
// An unknown token is always a hard error: the reader cannot skip a field of unknown size.
struct ExtraField
{
  ExtraFieldType type;
  std::string    name;
};
ExtraField parse_extra_field(const std::string& value);

// ---------------------------------------------------------------------------
// Header
// ---------------------------------------------------------------------------

struct QSMheader
{
  uint8_t  version_major = LIBQSM_VERSION_MAJOR;
  uint8_t  version_minor = LIBQSM_VERSION_MINOR;
  uint8_t  format        = 2;   // edge record layout: 0 = minimal, 1 = + arch, 2 = + distances
  std::string software;
  std::string created;
  std::string crs;
  uint64_t    n_nodes = 0;
  uint64_t    n_edges = 0;
  double      x_offset = 0.0;
  double      y_offset = 0.0;
  double      z_offset = 0.0;
  std::unordered_map<std::string, std::string> extra_keys;

  // Bounding box in global coordinates (mandatory; writer computes from nodes)
  double xmin = 0.0,  ymin = 0.0,  zmin = 0.0;
  double xmax = 0.0,  ymax = 0.0,  zmax = 0.0;

  // Named extra fields appended after the format-mandated bytes of each record
  std::vector<ExtraField> node_extra_fields;
  std::vector<ExtraField> edge_extra_fields;

  std::vector<std::string> messages;

  // Computed helpers: total extra bytes contributed by the declared fields.
  // Callers use spec.min_node_size + hdr.computed_node_extra_bytes(), etc.
  uint32_t computed_node_extra_bytes() const noexcept;
  uint32_t computed_edge_extra_bytes() const noexcept;
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
     << ", axis_id=" << e.axis_id
     << ", branch_order=" << static_cast<int>(e.branch_order)
     << ", subtree_length=" << e.subtree_length
     << ", distance_to_root=" << e.distance_to_root
     << ", quality=" << static_cast<int>(e.quality)
     << "}";
  return os;
}

// ---------------------------------------------------------------------------
// Format dispatch table
//
// Each format level owns four functions that speak *only* that format's
// binary layout.
//
//   read_*  reads exactly min_*_size mandatory bytes from the stream and
//           returns a fully-populated in-memory record; absent fields (those
//           not yet defined in that format) are left at their sentinel.
//
//   write_* writes exactly min_*_size bytes to the stream.
//
// The caller (read_binary / write_nodes / write_edges) is responsible for
// skipping any extra bytes beyond min_*_size that the file may carry (forward
// compatibility with a newer writer's extension fields).
// ---------------------------------------------------------------------------

struct FormatDispatch
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

struct QSMformatSpec
{
  uint8_t  format;
  uint32_t min_node_size;    // mandatory bytes per node record for this format
  uint32_t min_edge_size;    // mandatory bytes per edge record for this format
  FormatDispatch dispatch;   // reader/writer functors for this format
};

// ---------------------------------------------------------------------------
// Registry and resolver - defined in libqsm.cpp
//
// format_registry() returns the ordered table of known formats (lowest
// first).  Adding a new format means appending one row here; nothing else
// changes in the library core.
//
// resolve_format() selects the matching entry, or the highest known format
// when the file's FORMAT value is unrecognized (emitting a warning).
// ---------------------------------------------------------------------------

const std::vector<QSMformatSpec>& format_registry();
const QSMformatSpec& resolve_format(uint8_t format, std::string* warning_out = nullptr);

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
  uint8_t            get_format()        const noexcept { return header.format;        }
  const std::string& get_software()      const noexcept { return header.software;  }
  const std::string& get_created()       const noexcept { return header.created;   }
  const std::string& get_crs()           const noexcept { return header.crs;       }
  double             get_x_offset()      const noexcept { return header.x_offset;  }
  double             get_y_offset()      const noexcept { return header.y_offset;  }
  double             get_z_offset()      const noexcept { return header.z_offset;  }
  double             get_xmin()          const noexcept { return header.xmin; }
  double             get_ymin()          const noexcept { return header.ymin; }
  double             get_zmin()          const noexcept { return header.zmin; }
  double             get_xmax()          const noexcept { return header.xmax; }
  double             get_ymax()          const noexcept { return header.ymax; }
  double             get_zmax()          const noexcept { return header.zmax; }
  int                get_message_count() const noexcept { return static_cast<int>(header.messages.size()); }
  const std::string& get_message(int i)  const          { return header.messages.at(static_cast<std::size_t>(i)); }
  bool has_key(const std::string& key) const noexcept   { return header.extra_keys.find(key) != header.extra_keys.end(); }
  std::string get_key(const std::string& key) const noexcept { if (!has_key(key)) return {}; return header.extra_keys.at(key); }

  // Extra field declarations
  int               get_node_extra_count() const noexcept { return static_cast<int>(header.node_extra_fields.size()); }
  int               get_edge_extra_count() const noexcept { return static_cast<int>(header.edge_extra_fields.size()); }
  const ExtraField& get_node_extra_field(int i) const { return header.node_extra_fields.at(static_cast<std::size_t>(i)); }
  const ExtraField& get_edge_extra_field(int i) const { return header.edge_extra_fields.at(static_cast<std::size_t>(i)); }

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
  void set_format       (uint8_t f) noexcept             { header.format   = f;  }
  void set_software     (const std::string& s) noexcept  { header.software = s;  }
  void set_crs          (const std::string& c) noexcept  { header.crs      = c;  }
  void add_message      (const std::string& m) noexcept  { header.messages.push_back(m); }
  void set_origin(double x, double y, double z) noexcept { header.x_offset = x; header.y_offset = y; header.z_offset = z; }
  void add_key(const std::string& key, const std::string& value) noexcept { header.extra_keys[key] = value; }

  // Append an extra field declaration for node or edge records.
  // Fields are appended in call order, which is the binary storage order.
  // NOT supported in this implementation
  //void add_node_extra_field(ExtraFieldType type, const std::string& name) { header.node_extra_fields.push_back({type, name}); }
  //void add_edge_extra_field(ExtraFieldType type, const std::string& name) { header.edge_extra_fields.push_back({type, name}); }

  void add_node(const QSMnode& n) { nodes_.push_back(n); }
  void add_edge(const QSMedge& e) { edges_.push_back(e); }

  template <typename InputIt>
  void add_nodes(InputIt first, InputIt last) { nodes_.insert(nodes_.end(), first, last); }

  template <typename InputIt>
  void add_edges(InputIt first, InputIt last) { edges_.insert(edges_.end(), first, last); }

  void clear_records() noexcept { nodes_.clear(); edges_.clear(); }

  void write();

private:
  void write_header(std::ofstream& out, const QSMformatSpec& spec) const;
  void write_nodes (std::ofstream& out, const QSMformatSpec& spec) const;
  void write_edges (std::ofstream& out, const QSMformatSpec& spec) const;

  std::string filename_;

  QSMheader            header;
  std::vector<QSMnode> nodes_;
  std::vector<QSMedge> edges_;
};

} // namespace libqsm

#endif // LIBQSM_H

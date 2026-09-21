/**
 * @file libqsf.h
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
 *  QSF is a small, standalone, text-only manifest format that indexes a
 *  collection of .qsm files - conceptually similar to a GDAL VRT mosaic.
 *
 *  Two kinds are planned:
 *    - VIRTUAL  (implemented): the manifest only stores relative paths (and
 *      lightweight metadata such as id/name/bounding box) to external .qsm
 *      files sitting next to it. No tree data is duplicated.
 *    - EMBEDDED (reserved, not yet implemented): a future variant that will
 *      inline the full QSM data in the manifest itself, making it a fully
 *      self-contained "CSF" file. QSFreader already dispatches on KIND so
 *      this can be added later without breaking the format.
 *
 *  Write usage
 *
 *   libqsf::QSFwriter w("forest.qsf");
 *   w.set_software("Arbor");
 *   w.add_entry({"qsm/tree_1.qsm", 1, "tree_1", true, xmin,ymin,zmin,xmax,ymax,zmax});
 *   w.write();
 *
 *  Read usage
 *
 *   libqsf::QSFreader r("forest.qsf");
 *   for (const auto& e : r.entries()) { ... e.path, e.id, e.name ... }
 *
 */

#ifndef LIBQSF_H
#define LIBQSF_H

#include <cstdint>
#include <string>
#include <vector>

#define LIBQSF_VERSION_MAJOR 1
#define LIBQSF_VERSION_MINOR 0

#define LIBQSF_STRINGIFY2(x) #x
#define LIBQSF_STRINGIFY(x) LIBQSF_STRINGIFY2(x)

#define LIBQSF_VERSION_STRING              \
LIBQSF_STRINGIFY(LIBQSF_VERSION_MAJOR) "." \
LIBQSF_STRINGIFY(LIBQSF_VERSION_MINOR)

namespace libqsf {

// ---------------------------------------------------------------------------
// Kind of manifest
// ---------------------------------------------------------------------------

enum class QSFKind : uint8_t
{
  VIRTUAL  = 0,  // relative references to external .qsm files (implemented)
  EMBEDDED = 1   // full data inlined in the manifest (reserved, not yet implemented)
};

const char* kind_name(QSFKind k) noexcept;

// Throws std::runtime_error if 's' is not a recognized kind token.
QSFKind parse_kind(const std::string& s);

// ---------------------------------------------------------------------------
// One indexed QSM
// ---------------------------------------------------------------------------

struct QSFEntry
{
  std::string path;         // path to the .qsm file, relative to the .qsf file's own directory
  int32_t     id   = 0;
  std::string name;

  // Optional bounding box metadata (global coordinates), so callers can
  // inspect a forest's footprint without opening every referenced .qsm.
  bool   has_bbox = false;
  double xmin = 0.0, ymin = 0.0, zmin = 0.0;
  double xmax = 0.0, ymax = 0.0, zmax = 0.0;
};

// ---------------------------------------------------------------------------
// Header metadata
// ---------------------------------------------------------------------------

struct QSFheader
{
  uint8_t     version_major = LIBQSF_VERSION_MAJOR;
  uint8_t     version_minor = LIBQSF_VERSION_MINOR;
  QSFKind     kind          = QSFKind::VIRTUAL;
  std::string software;
  std::string created;
  std::string crs;
};

// ---------------------------------------------------------------------------
// QSFwriter
// ---------------------------------------------------------------------------

class QSFwriter
{
public:
  explicit QSFwriter(const std::string& filename);

  void set_version_major(uint8_t v) noexcept            { header.version_major = v; }
  void set_version_minor(uint8_t v) noexcept            { header.version_minor = v; }
  void set_kind          (QSFKind k) noexcept            { header.kind     = k;  }
  void set_software      (const std::string& s) noexcept { header.software = s;  }
  void set_crs           (const std::string& c) noexcept { header.crs      = c;  }
  void add_entry(const QSFEntry& e) { entries_.push_back(e); }

  void write() const;

private:
  std::string filename_;
  QSFheader header;
  std::vector<QSFEntry> entries_;
};

// ---------------------------------------------------------------------------
// QSFreader
// ---------------------------------------------------------------------------

class QSFreader
{
public:
  explicit QSFreader(const std::string& filename);

  uint8_t            get_version_major() const noexcept { return header.version_major; }
  uint8_t            get_version_minor() const noexcept { return header.version_minor; }
  QSFKind            get_kind()          const noexcept { return header.kind;          }
  const std::string& get_software()      const noexcept { return header.software;      }
  const std::string& get_created()       const noexcept { return header.created;       }
  const std::string& get_crs()           const noexcept { return header.crs;           }

  int                     entry_count() const noexcept { return static_cast<int>(entries_.size()); }
  const QSFEntry&         entry(int i)  const { return entries_.at(static_cast<std::size_t>(i)); }
  const std::vector<QSFEntry>& entries() const noexcept { return entries_; }

private:
  void parse();

  std::string filename_;
  QSFheader header;
  std::vector<QSFEntry> entries_;
};

} // namespace libqsf

#endif // LIBQSF_H

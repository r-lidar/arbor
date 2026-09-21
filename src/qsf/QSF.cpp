/**
 * @file QSF.cpp
 * Project: Arbor
 *
 * Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "QSF.h"
#include "libqsf.h"

#include <string>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <cctype>
#include <vector>
#include <array>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>

namespace fs = std::filesystem;

namespace arbor::qsm {

static std::array<double,3> compute_face_normal(const std::array<double,3>& v0, const std::array<double,3>& v1, const std::array<double,3>& v2)
{
  double ux = v1[0] - v0[0];
  double uy = v1[1] - v0[1];
  double uz = v1[2] - v0[2];

  double vx = v2[0] - v0[0];
  double vy = v2[1] - v0[1];
  double vz = v2[2] - v0[2];

  double nx = uy * vz - uz * vy;
  double ny = uz * vx - ux * vz;
  double nz = ux * vy - uy * vx;

  double len = std::sqrt(nx*nx + ny*ny + nz*nz);
  if (len <= 0.0) return {0.0, 0.0, 0.0};
  return { nx / len, ny / len, nz / len };
}

void QSF::add_qsm(const QSM& q)
{
  if (q.nodes().size() <= 1) return;
  qsm_[q.id] = q;
}

void QSF::write(const std::string& path, const std::string& format, bool binary) const
{
  if (path.empty())
    throw std::invalid_argument("QSF::write: output path is empty");

  if (format.empty())
    throw std::invalid_argument("QSF::write: format is empty");

  fs::path out_path(path);

  // Strip the leading '.' from the path's extension, if any.
  std::string ext = out_path.extension().string();
  if (!ext.empty() && ext.front() == '.')
    ext.erase(0, 1);

  // 'path' is treated as a single combined output FILE when it is not an
  // already-existing directory and its extension matches 'format'
  // (e.g. write("trees.obj", "obj")). Otherwise it is treated as a
  // directory that will receive one file per QSM (legacy behaviour).
  bool is_file_target = !ext.empty() && ext == format && !fs::is_directory(out_path);

  if (is_file_target)
  {
    if (format != "obj" && format != "ply" && format != "stl")
      throw std::invalid_argument("QSF::write: combined single-file output is only supported for 'obj', 'ply' and 'stl' formats, not '" + format + "'");

    fs::path parent = out_path.parent_path();
    if (!parent.empty() && !fs::exists(parent))
      fs::create_directories(parent);

    if (format == "obj")
      write_obj(out_path);
    else if (format == "ply")
      write_ply(out_path, binary);
    else
      write_stl(out_path, binary);

    return;
  }

  // Legacy behaviour: one file per QSM inside 'path/format/'
  fs::path base_dir(path);

  if (!fs::exists(base_dir))
  {
    fs::create_directories(base_dir);
  }
  else if (!fs::is_directory(base_dir))
  {
    throw std::runtime_error("QSF::write: path exists but is not a directory");
  }

  fs::path out_dir = base_dir / format;

  if (!fs::exists(out_dir))
    fs::create_directories(out_dir);

  std::unordered_map<int, fs::path> written_files;

  for (const auto& [key, qsm] : qsm_)
  {
    fs::path filename;
    if (qsm.name.empty())
      filename = out_dir / (std::to_string(qsm.id) + "." + format);
    else
      filename = out_dir / (qsm.name + "." + format);

    qsm.write(filename.string(), binary);
    written_files[key] = filename;
  }

  // A .qsm export additionally gets a sibling .qsf manifest (KIND VIRTUAL)
  // referencing every .qsm file just written, similar to a virtual raster
  // mosaic. Other formats (obj/ply/stl/csv/txt) are not indexed this way.
  if (format == "qsm")
    write_qsf_manifest(out_dir, written_files);
}

void QSF::write_obj(const fs::path& file) const
{
  std::ofstream out(file);
  if (!out.is_open()) throw std::runtime_error("Cannot open OBJ file.");

  out << std::fixed << std::setprecision(3);

  int vertex_offset = 0;

  for (const auto& [key, qsm] : qsm_)
  {
    std::vector<std::array<double,3>> vertices;
    std::vector<std::array<int,3>> faces;
    std::vector<int> node_ids;

    qsm.tmesh(vertices, faces, node_ids);

    // Each QSM is its own named object
    std::string label = qsm.name.empty() ? std::to_string(qsm.id) : qsm.name;
    for (auto& c : label)
      if (std::isspace(static_cast<unsigned char>(c))) c = '_';
    out << "o " << label << "\n";

    // Write vertices
    for (auto &v : vertices)
      out << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";

    // Write faces (OBJ uses 1-based indices)
    for (auto &f : faces)
      out << "f " << (f[0] + 1 + vertex_offset) << " " << (f[1] + 1 + vertex_offset) << " " << (f[2] + 1 + vertex_offset) << "\n";

    vertex_offset += static_cast<int>(vertices.size());
  }
}

void QSF::write_ply(const fs::path& file, bool binary) const
{
  struct Group
  {
    int id;
    std::string label;
    std::vector<std::array<double,3>> vertices;
    std::vector<std::array<int,3>> faces;
  };

  std::vector<Group> groups;
  groups.reserve(qsm_.size());

  size_t total_vertices = 0;
  size_t total_faces = 0;

  for (const auto& [key, qsm] : qsm_)
  {
    Group g;
    g.id = qsm.id;
    g.label = qsm.name.empty() ? std::to_string(qsm.id) : qsm.name;

    std::vector<int> node_ids;
    qsm.tmesh(g.vertices, g.faces, node_ids);

    total_vertices += g.vertices.size();
    total_faces += g.faces.size();
    groups.push_back(std::move(g));
  }

  if (binary)
  {
    std::ofstream out(file, std::ios::binary);
    if (!out.is_open()) throw std::runtime_error("Cannot open PLY file.");

    // header
    out << "ply\n";
    out << "format binary_little_endian 1.0\n";
    // one short comment per QSM keeps every group identifiable by id/name
    for (const auto& g : groups)
      out << "comment qsm " << g.id << " " << g.label << "\n";
    out << "element vertex " << total_vertices << "\n";
    out << "property double x\n";
    out << "property double y\n";
    out << "property double z\n";
    out << "element face " << total_faces << "\n";
    out << "property list uchar int vertex_indices\n";
    out << "property int qsm_id\n";
    out << "end_header\n";

    // write vertices
    for (const auto& g : groups)
    {
      for (auto &v : g.vertices)
      {
        out.write(reinterpret_cast<const char*>(&v[0]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&v[1]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&v[2]), sizeof(double));
      }
    }

    // write faces
    int vertex_offset = 0;
    for (const auto& g : groups)
    {
      for (auto &f : g.faces)
      {
        unsigned char nverts = 3; // triangle
        int a = f[0] + vertex_offset;
        int b = f[1] + vertex_offset;
        int c = f[2] + vertex_offset;
        out.write(reinterpret_cast<const char*>(&nverts), sizeof(unsigned char));
        out.write(reinterpret_cast<const char*>(&a), sizeof(int));
        out.write(reinterpret_cast<const char*>(&b), sizeof(int));
        out.write(reinterpret_cast<const char*>(&c), sizeof(int));
        out.write(reinterpret_cast<const char*>(&g.id), sizeof(int));
      }
      vertex_offset += static_cast<int>(g.vertices.size());
    }
  }
  else
  {
    std::ofstream out(file);
    if (!out.is_open()) throw std::runtime_error("Cannot open PLY file.");

    // header
    out << "ply\n";
    out << "format ascii 1.0\n";
    for (const auto& g : groups)
      out << "comment qsm " << g.id << " " << g.label << "\n";
    out << "element vertex " << total_vertices << "\n";
    out << "property double x\n";
    out << "property double y\n";
    out << "property double z\n";
    out << "element face " << total_faces << "\n";
    out << "property list uchar int vertex_indices\n";
    out << "property int qsm_id\n";
    out << "end_header\n";

    // write vertices
    for (const auto& g : groups)
      for (auto &v : g.vertices)
        out << std::fixed << std::setprecision(3) << v[0] << " " << v[1] << " " << v[2] << "\n";

    // write faces
    int vertex_offset = 0;
    for (const auto& g : groups)
    {
      for (auto &f : g.faces)
        out << "3 " << (f[0] + vertex_offset) << " " << (f[1] + vertex_offset) << " " << (f[2] + vertex_offset) << " " << g.id << "\n";
      vertex_offset += static_cast<int>(g.vertices.size());
    }
  }
}

void QSF::write_stl(const fs::path& file, bool binary) const
{
  struct Group
  {
    int id;
    std::string label;
    std::vector<std::array<double,3>> vertices;
    std::vector<std::array<int,3>> faces;
  };

  std::vector<Group> groups;
  groups.reserve(qsm_.size());

  for (const auto& [key, qsm] : qsm_)
  {
    Group g;
    g.id = qsm.id;
    g.label = qsm.name.empty() ? std::to_string(qsm.id) : qsm.name;

    std::vector<int> node_ids;
    qsm.tmesh(g.vertices, g.faces, node_ids);

    groups.push_back(std::move(g));
  }

  // STL is float-only and has no notion of geographic coordinates, so the
  // scene needs to be brought close to the origin. Unlike the single-QSM
  // writer, each tree can't be offset by its own root here: doing so would
  // collapse every tree onto (0,0,0) and make them all overlap once merged
  // into one file. Instead a single offset is computed once for the whole
  // scene - the centroid of every QSM's root node - and applied uniformly,
  // so trees keep their true positions relative to one another.
  double xoffset = 0.0, yoffset = 0.0, zoffset = 0.0;
  {
    double sx = 0.0, sy = 0.0, sz = 0.0;
    int n = 0;

    for (const auto& [key, qsm] : qsm_)
    {
      for (const auto& [nid, node] : qsm.nodes())
      {
        if (qsm.incoming_edges(nid).empty()) // root node
        {
          sx += node.x;
          sy += node.y;
          sz += node.z;
          ++n;
          break;
        }
      }
    }

    if (n > 0)
    {
      xoffset = sx / n;
      yoffset = sy / n;
      zoffset = sz / n;
    }
  }

  if (binary)
  {
    // Binary STL has no concept of multiple named solids or groups: it is
    // a flat header + triangle count + triangle list, with nothing to tag
    // a triangle as belonging to a particular QSM. All QSMs are therefore
    // merged here into a single, unnamed mesh. If keeping trees separate
    // and identifiable matters, use ASCII STL (binary = FALSE), or the
    // combined OBJ ('o' objects) / PLY ('qsm_id' property) writers instead.
    size_t total_faces = 0;
    for (const auto& g : groups)
      total_faces += g.faces.size();

    std::ofstream out(file, std::ios::binary);
    if (!out.is_open()) throw std::runtime_error("Cannot open STL file for writing: " + file.string());

    char header[80] = {};
    std::memcpy(header, "QSF binary STL (merged, unnamed)", 33);
    out.write(header, 80);

    uint32_t tri_count = static_cast<uint32_t>(total_faces);
    out.write(reinterpret_cast<const char*>(&tri_count), 4);

    for (const auto& g : groups)
    {
      for (const auto& f : g.faces)
      {
        const auto& v0 = g.vertices[f[0]];
        const auto& v1 = g.vertices[f[1]];
        const auto& v2 = g.vertices[f[2]];

        auto n = compute_face_normal(v0, v1, v2);
        float nf[3] = { float(n[0]), float(n[1]), float(n[2]) };
        out.write(reinterpret_cast<const char*>(nf), sizeof(nf));

        float fcoords[9] = {
          float(v0[0]-xoffset), float(v0[1]-yoffset), float(v0[2]-zoffset),
          float(v1[0]-xoffset), float(v1[1]-yoffset), float(v1[2]-zoffset),
          float(v2[0]-xoffset), float(v2[1]-yoffset), float(v2[2]-zoffset)
        };
        out.write(reinterpret_cast<const char*>(fcoords), sizeof(fcoords));

        uint16_t attr = 0;
        out.write(reinterpret_cast<const char*>(&attr), 2);
      }
    }
  }
  else
  {
    // ASCII STL: the format itself doesn't standardize multiple objects,
    // but the widely used convention is to concatenate several
    // "solid <name> ... endsolid <name>" blocks, one per object - mirroring
    // the 'o <name>' grouping used for OBJ. This keeps every QSM
    // distinguishable by name, but support for multiple solids in one file
    // is inconsistent across software: some readers import each block as a
    // separate named object, others only read the first block, and others
    // silently merge everything into one mesh. Check your target software
    // if this distinction matters.
    std::ofstream out(file);
    if (!out.is_open()) throw std::runtime_error("Cannot open STL file for writing: " + file.string());

    for (const auto& g : groups)
    {
      std::string label = g.label;
      for (auto& c : label)
        if (std::isspace(static_cast<unsigned char>(c))) c = '_';

      out << "solid " << label << "\n";
      for (const auto& f : g.faces)
      {
        const auto& v0 = g.vertices[f[0]];
        const auto& v1 = g.vertices[f[1]];
        const auto& v2 = g.vertices[f[2]];

        auto n = compute_face_normal(v0, v1, v2);
        out << "  facet normal " << n[0] << " " << n[1] << " " << n[2] << "\n";
        out << "    outer loop\n";
        out << "      vertex " << v0[0]-xoffset << " " << v0[1]-yoffset << " " << v0[2]-zoffset << "\n";
        out << "      vertex " << v1[0]-xoffset << " " << v1[1]-yoffset << " " << v1[2]-zoffset << "\n";
        out << "      vertex " << v2[0]-xoffset << " " << v2[1]-yoffset << " " << v2[2]-zoffset << "\n";
        out << "    endloop\n";
        out << "  endfacet\n";
      }
      out << "endsolid " << label << "\n";
    }
  }
}

void QSF::write_qsf_manifest(const fs::path& qsm_dir, const std::unordered_map<int, fs::path>& written_files) const
{
  // The manifest sits one level above the 'qsm/' directory it indexes, named
  // after that parent directory, e.g. writing to "forest" produces:
  //   forest/qsm/<name>.qsm  (one file per QSM, written by the caller)
  //   forest/forest.qsf      (this manifest, referencing the .qsm files above
  //                            by a path relative to its own directory)
  fs::path manifest_dir = qsm_dir.parent_path();
  std::string manifest_name = manifest_dir.filename().string();
  if (manifest_name.empty())
    manifest_name = "forest";
  fs::path manifest_file = manifest_dir / (manifest_name + ".qsf");

  libqsf::QSFwriter writer(manifest_file.string());
  writer.set_software("Arbor");
  writer.set_kind(libqsf::QSFKind::VIRTUAL);

  for (const auto& [key, qsm] : qsm_)
  {
    auto it = written_files.find(key);
    if (it == written_files.end()) continue; // should not happen

    libqsf::QSFEntry entry;
    entry.path = fs::relative(it->second, manifest_dir).string();
    entry.id   = qsm.id;
    entry.name = qsm.name;

    double xmin = std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double zmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymax = std::numeric_limits<double>::lowest();
    double zmax = std::numeric_limits<double>::lowest();

    for (const auto& [nid, node] : qsm.nodes())
    {
      xmin = std::min(xmin, node.x); xmax = std::max(xmax, node.x);
      ymin = std::min(ymin, node.y); ymax = std::max(ymax, node.y);
      zmin = std::min(zmin, node.z); zmax = std::max(zmax, node.z);
    }

    if (!qsm.nodes().empty())
    {
      entry.has_bbox = true;
      entry.xmin = xmin; entry.ymin = ymin; entry.zmin = zmin;
      entry.xmax = xmax; entry.ymax = ymax; entry.zmax = zmax;
    }

    writer.add_entry(entry);
  }

  writer.write();
}

QSF QSF::read(const std::string& path)
{
  fs::path manifest_file(path);

  if (!fs::exists(manifest_file))
    throw std::runtime_error("QSF::read: file not found: " + path);

  libqsf::QSFreader reader(manifest_file.string());

  // KIND dispatch. libqsf::QSFreader already refuses to parse an EMBEDDED
  // manifest's entries, but the check is repeated here (defensively, and to
  // keep the branch structure explicit and easy to extend once an EMBEDDED
  // reader is implemented).
  QSF result;

  if (reader.get_kind() == libqsf::QSFKind::VIRTUAL)
  {
    fs::path manifest_dir = manifest_file.parent_path();

    for (const auto& entry : reader.entries())
    {
      fs::path qsm_file = manifest_dir / entry.path;

      if (!fs::exists(qsm_file))
        throw std::runtime_error("QSF::read: referenced .qsm file not found: " + qsm_file.string() + " (from manifest " + path + ")");

      QSM qsm;
      qsm.read(qsm_file.string());

      // The manifest's id/name take precedence over whatever the .qsm file
      // itself stores, since the manifest is the authoritative index.
      qsm.id = entry.id;
      if (!entry.name.empty())
        qsm.name = entry.name;

      result.add_qsm(qsm);
    }
  }
  else if (reader.get_kind() == libqsf::QSFKind::EMBEDDED)
  {
    // Reserved for a future fully self-contained QSF variant. libqsf already
    // throws while parsing an EMBEDDED manifest's body, so this branch is
    // effectively unreachable today; it is kept to make the intended
    // extension point explicit.
    throw std::runtime_error("QSF::read: KIND EMBEDDED is not yet supported.");
  }
  else
  {
    throw std::runtime_error("QSF::read: unsupported manifest KIND.");
  }

  return result;
}

}
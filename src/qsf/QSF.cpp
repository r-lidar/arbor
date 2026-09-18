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

#include <array>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <string>
#include <filesystem>
#include <stdexcept>
#include <vector>

namespace fs = std::filesystem;

namespace arbor::qsm {

namespace {

std::string to_lower(std::string s)
{
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
  return s;
}

} // anonymous namespace

void QSF::add_qsm(const QSM& q)
{
  if (q.nodes().size() <= 1) return;
  qsm_[q.id] = q;
}

void QSF::write_single_file(const std::string& filename, const std::string& format, bool binary) const
{
  if (format != "obj" && format != "ply")
    throw std::runtime_error("QSF::write: writing all QSM into a single file is only supported for .obj and .ply formats");

  // OBJ has no official binary variant: it is always a plain-text format,
  // exactly as for a single QSM (see QSM::write_obj). The 'binary' flag is
  // only meaningful for PLY.

  // Merge the meshes of every QSM into a single set of vertices/faces, but
  // keep track of which QSM each vertex/face came from so that each QSM can
  // be written back out as its own named object (OBJ "o" group) / tagged
  // face (PLY "object_id" property), instead of a single flattened blob.
  std::vector<std::array<double,3>> vertices;
  std::vector<std::array<int,3>> faces;
  std::vector<int> face_object_id;
  std::vector<std::string> object_names;
  std::vector<std::pair<size_t,size_t>> object_vertex_range; // [start, end)
  std::vector<std::pair<size_t,size_t>> object_face_range;   // [start, end)

  for (const auto& [key, qsm] : qsm_)
  {
    std::vector<std::array<double,3>> v;
    std::vector<std::array<int,3>> f;
    std::vector<int> node_ids;

    qsm.tmesh(v, f, node_ids);

    int voffset = static_cast<int>(vertices.size());
    size_t vstart = vertices.size();
    size_t fstart = faces.size();

    vertices.insert(vertices.end(), v.begin(), v.end());

    int object_id = static_cast<int>(object_names.size());
    for (const auto& face : f)
    {
      faces.push_back({ face[0] + voffset, face[1] + voffset, face[2] + voffset });
      face_object_id.push_back(object_id);
    }

    object_names.push_back(qsm.name.empty() ? std::to_string(qsm.id) : qsm.name);
    object_vertex_range.emplace_back(vstart, vertices.size());
    object_face_range.emplace_back(fstart, faces.size());
  }

  if (format == "obj")
  {
    std::ofstream out(filename);
    if (!out.is_open()) throw std::runtime_error("Cannot open OBJ file: " + filename);

    out << std::fixed << std::setprecision(3);

    for (size_t o = 0; o < object_names.size(); ++o)
    {
      out << "o " << object_names[o] << "\n";

      auto [vstart, vend] = object_vertex_range[o];
      for (size_t i = vstart; i < vend; ++i)
        out << "v " << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << "\n";

      auto [fstart, fend] = object_face_range[o];
      for (size_t i = fstart; i < fend; ++i)
        out << "f " << (faces[i][0]+1) << " " << (faces[i][1]+1) << " " << (faces[i][2]+1) << "\n";
    }
  }
  else // ply
  {
    if (binary)
    {
      std::ofstream out(filename, std::ios::binary);
      if (!out.is_open()) throw std::runtime_error("Cannot open PLY file: " + filename);

      out << "ply\n";
      out << "format binary_little_endian 1.0\n";
      for (size_t o = 0; o < object_names.size(); ++o)
        out << "comment object " << o << " " << object_names[o] << "\n";
      out << "element vertex " << vertices.size() << "\n";
      out << "property double x\n";
      out << "property double y\n";
      out << "property double z\n";
      out << "element face " << faces.size() << "\n";
      out << "property list uchar int vertex_indices\n";
      out << "property int object_id\n";
      out << "end_header\n";

      for (const auto& v : vertices)
      {
        out.write(reinterpret_cast<const char*>(&v[0]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&v[1]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&v[2]), sizeof(double));
      }

      for (size_t i = 0; i < faces.size(); ++i)
      {
        const auto& f = faces[i];
        unsigned char nverts = 3;
        out.write(reinterpret_cast<const char*>(&nverts), sizeof(unsigned char));
        out.write(reinterpret_cast<const char*>(&f[0]), sizeof(int));
        out.write(reinterpret_cast<const char*>(&f[1]), sizeof(int));
        out.write(reinterpret_cast<const char*>(&f[2]), sizeof(int));
        out.write(reinterpret_cast<const char*>(&face_object_id[i]), sizeof(int));
      }
    }
    else
    {
      std::ofstream out(filename);
      if (!out.is_open()) throw std::runtime_error("Cannot open PLY file: " + filename);

      out << "ply\n";
      out << "format ascii 1.0\n";
      for (size_t o = 0; o < object_names.size(); ++o)
        out << "comment object " << o << " " << object_names[o] << "\n";
      out << "element vertex " << vertices.size() << "\n";
      out << "property double x\n";
      out << "property double y\n";
      out << "property double z\n";
      out << "element face " << faces.size() << "\n";
      out << "property list uchar int vertex_indices\n";
      out << "property int object_id\n";
      out << "end_header\n";

      for (const auto& v : vertices)
        out << std::fixed << std::setprecision(3) << v[0] << " " << v[1] << " " << v[2] << "\n";

      for (size_t i = 0; i < faces.size(); ++i)
      {
        const auto& f = faces[i];
        out << "3 " << f[0] << " " << f[1] << " " << f[2] << " " << face_object_id[i] << "\n";
      }
    }
  }
}

void QSF::write(const std::string& dir, const std::string& format, bool binary) const
{
  if (dir.empty())
    throw std::invalid_argument("QSF::write: output directory is empty");

  if (format.empty())
    throw std::invalid_argument("QSF::write: format is empty");

  std::string fmt = to_lower(format);

  // If the output path itself carries an extension matching the requested
  // format, treat it as a single output file gathering every QSM instead of
  // a directory containing one file per QSM.
  fs::path out_path(dir);
  if (out_path.has_extension())
  {
    std::string ext = to_lower(out_path.extension().string());
    if (!ext.empty() && ext.front() == '.') ext.erase(ext.begin());

    if (ext == fmt)
    {
      fs::path parent = out_path.parent_path();
      if (!parent.empty() && !fs::exists(parent))
        fs::create_directories(parent);

      write_single_file(out_path.string(), fmt, binary);
      return;
    }
  }

  fs::path base_dir(dir);

  // Create base directory if needed
  if (!fs::exists(base_dir))
  {
    fs::create_directories(base_dir);
  }
  else if (!fs::is_directory(base_dir))
  {
    throw std::runtime_error("QSF::write: path exists but is not a directory");
  }

  // Optional subfolder per format
  fs::path out_dir = base_dir / format;

  if (!fs::exists(out_dir))
    fs::create_directories(out_dir);

  // Write each QSM
  for (const auto& [key, qsm] : qsm_)
  {
    fs::path filename;
    if (qsm.name.empty())
      filename = out_dir / (std::to_string(qsm.id) + "." + format);
    else
      filename = out_dir / (qsm.name + "." + format);

    qsm.write(filename.string(), binary);
  }
}

}

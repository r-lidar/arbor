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

  // Merge the meshes of every QSM into a single set of vertices/faces
  std::vector<std::array<double,3>> vertices;
  std::vector<std::array<int,3>> faces;

  for (const auto& [key, qsm] : qsm_)
  {
    std::vector<std::array<double,3>> v;
    std::vector<std::array<int,3>> f;
    std::vector<int> node_ids;

    qsm.tmesh(v, f, node_ids);

    int offset = static_cast<int>(vertices.size());
    vertices.insert(vertices.end(), v.begin(), v.end());
    for (const auto& face : f)
      faces.push_back({ face[0] + offset, face[1] + offset, face[2] + offset });
  }

  if (format == "obj")
  {
    std::ofstream out(filename);
    if (!out.is_open()) throw std::runtime_error("Cannot open OBJ file: " + filename);

    for (const auto& v : vertices)
      out << "v " << std::fixed << std::setprecision(3) << v[0] << " " << v[1] << " " << v[2] << "\n";

    for (const auto& f : faces)
      out << "f " << (f[0]+1) << " " << (f[1]+1) << " " << (f[2]+1) << "\n";
  }
  else // ply
  {
    if (binary)
    {
      std::ofstream out(filename, std::ios::binary);
      if (!out.is_open()) throw std::runtime_error("Cannot open PLY file: " + filename);

      out << "ply\n";
      out << "format binary_little_endian 1.0\n";
      out << "element vertex " << vertices.size() << "\n";
      out << "property double x\n";
      out << "property double y\n";
      out << "property double z\n";
      out << "element face " << faces.size() << "\n";
      out << "property list uchar int vertex_indices\n";
      out << "end_header\n";

      for (const auto& v : vertices)
      {
        out.write(reinterpret_cast<const char*>(&v[0]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&v[1]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&v[2]), sizeof(double));
      }

      for (const auto& f : faces)
      {
        unsigned char nverts = 3;
        out.write(reinterpret_cast<const char*>(&nverts), sizeof(unsigned char));
        out.write(reinterpret_cast<const char*>(&f[0]), sizeof(int));
        out.write(reinterpret_cast<const char*>(&f[1]), sizeof(int));
        out.write(reinterpret_cast<const char*>(&f[2]), sizeof(int));
      }
    }
    else
    {
      std::ofstream out(filename);
      if (!out.is_open()) throw std::runtime_error("Cannot open PLY file: " + filename);

      out << "ply\n";
      out << "format ascii 1.0\n";
      out << "element vertex " << vertices.size() << "\n";
      out << "property double x\n";
      out << "property double y\n";
      out << "property double z\n";
      out << "element face " << faces.size() << "\n";
      out << "property list uchar int vertex_indices\n";
      out << "end_header\n";

      for (const auto& v : vertices)
        out << std::fixed << std::setprecision(3) << v[0] << " " << v[1] << " " << v[2] << "\n";

      for (const auto& f : faces)
        out << "3 " << f[0] << " " << f[1] << " " << f[2] << "\n";
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

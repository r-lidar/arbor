/**
 * @file QSM_io.cpp
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

#include "QSM.h"
#include "libqsm.h"

#include <fstream>
#include <iomanip> // std::setprecision
#include <algorithm>
#include <cmath>
#include <cstring> // std::memset
#include <cstdint>
#include <ctime>   // std::time, std::strftime, std::gmtime
#include <array>

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

void QSM::write(const std::string& filename, bool binary) const
{
  if (edges().size() == 0) return;

  // Find extension
  auto pos = filename.find_last_of('.');
  if (pos == std::string::npos) throw std::runtime_error("Filename has no extension: " + filename);

  std::string ext = filename.substr(pos + 1);

  // Normalize extension
  std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return std::tolower(c); });

  if (ext == "qsm")
    write_qsm(filename);
  else if (ext == "obj")
    write_obj(filename);
  else if (ext == "stl")
    write_stl(filename, binary);
  else if (ext == "ply")
    write_ply(filename, binary);
  else if (ext == "csv" || ext == "txt")
    write_csv(filename);
  else
    throw std::runtime_error("Unknown file extension: " + ext + ". Supported: .qbf, .obj, .ply, .stl, .csv, .txt");
}

void QSM::write_qsm(const std::string& filename) const
{
  libqsm::QSMwriter writer(filename);
  writer.set_software("Arbor");
  writer.set_version_minor(2);

  double xoffset = 0.0, yoffset = 0.0, zoffset = 0.0;
  {
    const NodeID root = find_root_node();
    if (root != -1)
    {
      const QSMNode& r = node(root);
      xoffset = r.x;
      yoffset = r.y;
      zoffset = r.z;
    }
  }
  writer.set_origin(xoffset, yoffset, zoffset);
  for (const auto& elem : nodes()) writer.add_node(elem.second);
  for (const auto& elem : edges()) writer.add_edge(elem.second.data);
  writer.write();
}


void QSM::write_ply(const std::string& filename, bool binary) const
{
  std::vector<std::array<double,3>> vertices;
  std::vector<std::array<int,3>> faces;
  std::vector<int> node_ids;

  tmesh(vertices, faces, node_ids);

  if (binary)
  {
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) throw std::runtime_error("Cannot open PLY file.");

    // header
    out << "ply\n";
    out << "format binary_little_endian 1.0\n";
    out << "element vertex " << vertices.size() << "\n";
    out << "property double x\n";
    out << "property double y\n";
    out << "property double z\n";
    out << "element face " << faces.size() << "\n";
    out << "property list uchar int vertex_indices\n";
    out << "end_header\n";

    // write vertices
    for (auto &v : vertices)
    {
      out.write(reinterpret_cast<const char*>(&v[0]), sizeof(double));
      out.write(reinterpret_cast<const char*>(&v[1]), sizeof(double));
      out.write(reinterpret_cast<const char*>(&v[2]), sizeof(double));
    }

    // write faces
    for (auto &f : faces)
    {
      unsigned char nverts = 3; // triangle
      out.write(reinterpret_cast<const char*>(&nverts), sizeof(unsigned char));
      out.write(reinterpret_cast<const char*>(&f[0]), sizeof(int));
      out.write(reinterpret_cast<const char*>(&f[1]), sizeof(int));
      out.write(reinterpret_cast<const char*>(&f[2]), sizeof(int));
    }
  }
  else
  {
    std::ofstream out(filename);
    if (!out.is_open()) throw std::runtime_error("Cannot open PLY file.");

    // header
    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << vertices.size() << "\n";
    out << "property double x\n";
    out << "property double y\n";
    out << "property double z\n";
    out << "element face " << faces.size() << "\n";
    out << "property list uchar int vertex_indices\n";
    out << "end_header\n";

    // write vertices
    for (auto &v : vertices)
      out << std::fixed << std::setprecision(3) << v[0] << " " << v[1] << " " << v[2] << "\n";

    // write faces
    for (auto &f : faces)
      out << "3 " << f[0] << " " << f[1] << " " << f[2] << "\n";
  }
}

void QSM::write_obj(const std::string& filename) const
{
  std::vector<std::array<double,3>> vertices;
  std::vector<std::array<int,3>> faces;
  std::vector<int> node_ids;

  tmesh(vertices, faces, node_ids);

  std::ofstream out(filename);
  if (!out.is_open()) throw std::runtime_error("Cannot open OBJ file.");

  // Write vertices
  for (auto &v : vertices)
    out << "v " << std::fixed << std::setprecision(3) << v[0] << " " << v[1] << " " << v[2] << "\n";

  // Write faces (OBJ uses 1-based indices)
  for (auto &f : faces)
    out << "f " << (f[0]+1) << " " << (f[1]+1) << " " << (f[2]+1) << "\n";
}

void QSM::write_stl(const std::string& filename, bool binary) const
{
  std::vector<std::array<double,3>> vertices;
  std::vector<std::array<int,3>> faces;
  std::vector<int> node_ids;

  tmesh(vertices, faces, node_ids);

  // Determine offset from the root node
  // because STL is float only and does not support geographic coordinates
  double xoffset = 0.0, yoffset = 0.0, zoffset = 0.0;

  for (const auto& [nid, n] : nodes())
  {
    if (incoming_edges(nid).empty()) // root node
    {
      xoffset = n.x;
      yoffset = n.y;
      zoffset = n.z;
      break;
    }
  }

  if (binary)
  {
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open())
      throw std::runtime_error("Cannot open STL file for writing: " + filename);

    char header[80] = {};
    std::memcpy(header, "QSM binary STL", 14);
    out.write(header, 80);

    // Number of triangles
    uint32_t tri_count = static_cast<uint32_t>(faces.size());
    out.write(reinterpret_cast<const char*>(&tri_count), 4);

    // Write triangles
    for (const auto& f : faces)
    {
      const auto& v0 = vertices[f[0]];
      const auto& v1 = vertices[f[1]];
      const auto& v2 = vertices[f[2]];

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
  else
  {
    // ASCII STL
    std::ofstream out(filename);
    if (!out.is_open())
      throw std::runtime_error("Cannot open STL file for writing: " + filename);

    out << "solid QSM\n";
    for (const auto& f : faces)
    {
      const auto& v0 = vertices[f[0]];
      const auto& v1 = vertices[f[1]];
      const auto& v2 = vertices[f[2]];

      auto n = compute_face_normal(v0, v1, v2);
      out << "  facet normal " << n[0] << " " << n[1] << " " << n[2] << "\n";
      out << "    outer loop\n";
      out << "      vertex " << v0[0]-xoffset << " " << v0[1]-yoffset << " " << v0[2]-zoffset << "\n";
      out << "      vertex " << v1[0]-xoffset << " " << v1[1]-yoffset << " " << v1[2]-zoffset << "\n";
      out << "      vertex " << v2[0]-xoffset << " " << v2[1]-yoffset << " " << v2[2]-zoffset << "\n";
      out << "    endloop\n";
      out << "  endfacet\n";
    }
    out << "endsolid QSM\n";
  }
}


void QSM::write_csv(const std::string& filename) const
{
  std::ofstream out(filename);
  if (!out.is_open())
    throw std::runtime_error("Cannot open CSV file: " + filename);

  // Header
  out << "startX,startY,startZ,endX,endY,endZ,id,source,axis_id,branch_order,radius,length,volume,subtree_length";
  out << std::endl;

  out << std::fixed << std::setprecision(3);

  // ---- SORT EDGES BY id ----
  std::vector<const EdgeInfo*> sorted;
  sorted.reserve(edges().size());

  for (const auto& kv : edges())
    sorted.push_back(&kv.second);

  std::sort(sorted.begin(), sorted.end(),
            [](const EdgeInfo* a, const EdgeInfo* b)
            {
              return a->data.id < b->data.id;
            });

  // ---- WRITE DATA IN ORDER ----
  for (const EdgeInfo* e : sorted)
  {
    const QSMNode& src = node(e->source);
    const QSMNode& tgt = node(e->target);
    const QSMEdge& c   = e->data;

    out << src.x           << ","
        << src.y           << ","
        << src.z           << ","
        << tgt.x           << ","
        << tgt.y           << ","
        << tgt.z           << ","
        << c.id            << ","
        << c.source        << ","
        << c.axis_id       << ","
        << c.branch_order  << ","
        << c.radius        << ","
        << c.length(src,tgt) << ","
        << std::fixed << std::setprecision(5)
        << c.volume(src,tgt) << ","
        << std::fixed << std::setprecision(3)
        << c.subtree_length
        << std::endl;
  }
}

} // namespace arbor::qsm

#include "QSM.h"
#include <fstream>
#include <iomanip> // std::setprecision
#include <algorithm>

void QSM::write(const std::string& filename, int resolution) const
{
  // Find extension
  auto pos = filename.find_last_of('.');
  if (pos == std::string::npos) throw std::runtime_error("Filename has no extension: " + filename);

  std::string ext = filename.substr(pos + 1);

  // Normalize extension
  std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return std::tolower(c); });

  if (ext == "obj")
    write_obj(filename, resolution);
  else if (ext == "ply")
    write_ply(filename, resolution);
  else if (ext == "csv" || ext == "txt")
    write_csv(filename);
  else
    throw std::runtime_error("Unknown file extension: " + ext + ". Supported: obj, ply, csv, txt");
}

void QSM::write_ply(const std::string& filename, int resolution) const
{
  std::vector<std::array<double,3>> vertices;
  std::vector<std::array<int,3>> faces;

  build_mesh(vertices, faces, resolution);

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

void QSM::write_obj(const std::string& filename, int resolution) const
{
  std::vector<std::array<double,3>> vertices;
  std::vector<std::array<int,3>> faces;

  build_mesh(vertices, faces, resolution);

  std::ofstream out(filename);
  if (!out.is_open()) throw std::runtime_error("Cannot open OBJ file.");

  // Write vertices
  for (auto &v : vertices)
    out << "v " << std::fixed << std::setprecision(3) << v[0] << " " << v[1] << " " << v[2] << "\n";

  // Write faces (OBJ uses 1-based indices)
  for (auto &f : faces)
    out << "f " << (f[0]+1) << " " << (f[1]+1) << " " << (f[2]+1) << "\n";
}

void QSM::write_csv(const std::string& filename) const
{
  std::ofstream out(filename);
  if (!out.is_open())
    throw std::runtime_error("Cannot open CSV file: " + filename);

  // Header
  out <<   "startX,startY,startZ,endX,endY,endZ,cyl_ID,parent_ID,axis_ID,branch_order,radius,length,volume,subtree_length";
  out << std::endl;

  out << std::fixed << std::setprecision(3);

  // ---- SORT CYLINDERS BY cyl_ID ----
  std::vector<const QSMcylinder*> sorted;
  sorted.reserve(cylinders_.size());
  for (const auto& kv : cylinders_) sorted.push_back(&kv.second);
  std::sort(sorted.begin(), sorted.end(), [](const QSMcylinder* a, const QSMcylinder* b) {return a->cyl_ID < b->cyl_ID; });

  // ---- WRITE DATA IN ORDER ----
  for (auto* c : sorted)
  {
    out << c->startX              << ","
        << c->startY              << ","
        << c->startZ              << ","
        << c->endX                << ","
        << c->endY                << ","
        << c->endZ                << ","
        << c->cyl_ID              << ","
        << c->parent_ID           << ","
        << c->axis_ID             << ","
        << c->branch_order        << ","
        << c->radius              << ","
        << c->length()            << ","
        << c->volume()            << ","
        << c->subtree_length
        << std::endl;
  }
}

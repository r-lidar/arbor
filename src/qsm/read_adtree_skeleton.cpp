#include <Rcpp.h>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>

using namespace Rcpp;
using namespace std;

struct Vertex
{
  float x, y, z, radius;
};

struct Edge
{
  uint32_t v1, v2;
};

DataFrame read_adtree_skeleton(std::string filename)
{
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open()) stop("Cannot open file");

  std::string line;
  size_t num_vertices = 0, num_edges = 0;

  // === 1. Read Header ===
  while (std::getline(file, line)) {
    size_t pos;
    if ((pos = line.find("element vertex")) != std::string::npos) {
      num_vertices = std::stoul(line.substr(pos + 14));
    } else if ((pos = line.find("element edge")) != std::string::npos) {
      num_edges = std::stoul(line.substr(pos + 13));
    } else if (line == "end_header") {
      break;
    }
  }

  // === 2. Read Vertices ===
  std::vector<Vertex> vertices(num_vertices);
  for (size_t i = 0; i < num_vertices; ++i) {
    file.read(reinterpret_cast<char*>(&vertices[i].x), sizeof(float));
    file.read(reinterpret_cast<char*>(&vertices[i].y), sizeof(float));
    file.read(reinterpret_cast<char*>(&vertices[i].z), sizeof(float));
    file.read(reinterpret_cast<char*>(&vertices[i].radius), sizeof(float));
  }

  // === 3. Read Edges ===
  std::vector<Edge> edges;

  for (size_t i = 0; i < num_edges; ++i) {
    uint32_t list_size;
    file.read(reinterpret_cast<char*>(&list_size), sizeof(uint32_t));

    if (file.fail()) stop("Error reading edge list size");

    if (list_size < 2) continue;

    std::vector<int32_t> indices(list_size);
    file.read(reinterpret_cast<char*>(indices.data()), list_size * sizeof(int32_t));

    if (file.fail()) stop("Error reading edge vertex indices");

    if (list_size == 2) {
      edges.push_back({static_cast<uint32_t>(indices[0]), static_cast<uint32_t>(indices[1])});
    } else if (list_size > 2) {
      for (size_t j = 0; j < list_size - 1; ++j) {
        edges.push_back({static_cast<uint32_t>(indices[j]), static_cast<uint32_t>(indices[j + 1])});
      }
    } else {
      Rcout << "Skipping edge with less than 2 vertices" << std::endl;
    }
  }

  // === 4. Build R DataFrame ===
  std::vector<float> x1, y1, z1, r1, x2, y2, z2, r2, id, parentID;
  x1.reserve(num_edges); y1.reserve(num_edges); z1.reserve(num_edges); r1.reserve(num_edges);
  x2.reserve(num_edges); y2.reserve(num_edges); z2.reserve(num_edges); r2.reserve(num_edges);
  id.reserve(num_edges); parentID.reserve(num_edges);

  for (const auto& edge : edges) {
    if (edge.v1 >= vertices.size() || edge.v2 >= vertices.size()) {
      Rcout << "Invalid vertex indices: " << edge.v1 << ", " << edge.v2 << std::endl;
      continue; // Skip this edge
    }

    const Vertex& v1 = vertices[edge.v1];
    const Vertex& v2 = vertices[edge.v2];

    x1.push_back(v1.x); y1.push_back(v1.y); z1.push_back(v1.z); r1.push_back(v1.radius);
    x2.push_back(v2.x); y2.push_back(v2.y); z2.push_back(v2.z); r2.push_back(v2.radius);
    id.push_back(edge.v2); parentID.push_back(edge.v1);
  }

  return DataFrame::create(
    _["startX"] = x1, _["startY"] = y1, _["startZ"] = z1, _["radius1"] = r1,
      _["endX"] = x2, _["endY"] = y2, _["endZ"] = z2, _["radius2"] = r2,
        _["cyl_ID"] = id, _["parent_ID"] = parentID
  );
}

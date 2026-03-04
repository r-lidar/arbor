#ifndef GRID3D_H
#define GRID3D_H

#include <vector>
#include <unordered_map>
#include <cstdint>

#include "PointCloud.h"

struct Offset
{
  int dx, dy, dz;
};

class Grid3D
{
public:
  Grid3D(const PointCloud& pc, double res);

  // Grid dimensions
  int64_t ncols, nrows, nlayers;

  // Bounding box
  double xmin, ymin, zmin;
  double xmax, ymax, zmax;

  // Inverse resolution
  double inv_xres, inv_yres, inv_zres;

  int npoints;

  std::vector<int> connected_components(int connectivity);

private:
  const PointCloud& pc;

  int64_t get_grid_index(int64_t r, int64_t c, int64_t l) const;
  void decode_grid_index(int64_t idx, int64_t& r, int64_t& c, int64_t& l) const;
  std::vector<Offset> get_neighbor_offsets(int connectivity) const;
};

#endif // GRID3D_H

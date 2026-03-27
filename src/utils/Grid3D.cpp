#include "Grid3D.h"
#include "services.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>

Grid3D::Grid3D(const PointCloud& pc, double res) : pc(pc)
{
  npoints = pc.size();

  if (npoints <= 0)
  {
    xmin = xmax = ymin = ymax = zmin = zmax = 0.0;
    ncols = nrows = nlayers = 0;
    inv_xres = inv_yres = inv_zres = 0.0;
    return;
  }

  xmin = xmax = pc.get_x(0);
  ymin = ymax = pc.get_y(0);
  zmin = zmax = pc.get_z(0);

  for (int i = 1; i < npoints; ++i)
  {
    xmin = std::min(xmin, pc.get_x(i));
    xmax = std::max(xmax, pc.get_x(i));
    ymin = std::min(ymin, pc.get_y(i));
    ymax = std::max(ymax, pc.get_y(i));
    zmin = std::min(zmin, pc.get_z(i));
    zmax = std::max(zmax, pc.get_z(i));
  }

  double inv_res = 1.0 / res;
  inv_xres = inv_yres = inv_zres = inv_res;

  xmin = std::floor(xmin * inv_xres) * res;
  ymin = std::floor(ymin * inv_yres) * res;
  zmin = std::floor(zmin * inv_zres) * res;

  ncols   = static_cast<int64_t>((xmax - xmin) * inv_xres) + 1;
  nrows   = static_cast<int64_t>((ymax - ymin) * inv_yres) + 1;
  nlayers = static_cast<int64_t>((zmax - zmin) * inv_zres) + 1;
}

int64_t Grid3D::get_grid_index(int64_t r, int64_t c, int64_t l) const
{
  return l * (nrows * ncols) + r * ncols + c;
}

void Grid3D::decode_grid_index(int64_t idx, int64_t& r, int64_t& c, int64_t& l) const
{
  int64_t area = nrows * ncols;
  l = idx / area;
  int64_t rem = idx % area;
  r = rem / ncols;
  c = rem % ncols;
}

std::vector<Offset> Grid3D::get_neighbor_offsets(int connectivity) const
{
  std::vector<Offset> offsets;

  for (int dz = -1; dz <= 1; ++dz)
  {
    for (int dy = -1; dy <= 1; ++dy)
    {
      for (int dx = -1; dx <= 1; ++dx)
      {
        if (dx == 0 && dy == 0 && dz == 0)
          continue;

        int sum_abs = std::abs(dx) + std::abs(dy) + std::abs(dz);

        bool ok = false;
        if (connectivity == 6  && sum_abs == 1) ok = true;
        if (connectivity == 18 && sum_abs <= 2) ok = true;
        if (connectivity == 26) ok = true;

        if (ok)
          offsets.push_back({dx, dy, dz});
      }
    }
 }

 return offsets;
}

std::vector<int> Grid3D::connected_components(int connectivity)
{
  if (npoints == 0)
    return {};

  std::vector<int64_t> point_to_voxel(npoints);
  std::unordered_map<int64_t, int> voxel_labels;
  voxel_labels.reserve(npoints / 2);

  for (int i = 0; i < npoints; ++i)
  {
    int64_t c = static_cast<int64_t>((pc.get_x(i) - xmin) * inv_xres);
    int64_t r = static_cast<int64_t>((pc.get_y(i) - ymin) * inv_yres);
    int64_t l = static_cast<int64_t>((pc.get_z(i) - zmin) * inv_zres);

    c = std::clamp(c, int64_t(0), ncols   - 1);
    r = std::clamp(r, int64_t(0), nrows   - 1);
    l = std::clamp(l, int64_t(0), nlayers - 1);

    int64_t idx = get_grid_index(r, c, l);
    point_to_voxel[i] = idx;
    voxel_labels[idx] = 0;
  }

  auto neighbors = get_neighbor_offsets(connectivity);

  int current_label = 0;
  std::vector<int64_t> stack;
  stack.reserve(1024);

  auto prog = ServiceLocator::make_progress(voxel_labels.size(), "Connected components", 0.25);

  for (auto& kv : voxel_labels)
  {
    if (kv.second != 0)
      continue;

    ++current_label;
    kv.second = current_label;
    stack.push_back(kv.first);
    prog->tick();

    while (!stack.empty())
    {
      int64_t curr = stack.back();
      stack.pop_back();

      int64_t r, c, l;
      decode_grid_index(curr, r, c, l);

      for (const auto& off : neighbors)
      {
        int64_t nr = r + off.dy;
        int64_t nc = c + off.dx;
        int64_t nl = l + off.dz;

        if (nr < 0 || nr >= nrows ||
            nc < 0 || nc >= ncols ||
            nl < 0 || nl >= nlayers)
          continue;

        int64_t nidx = get_grid_index(nr, nc, nl);
        auto it = voxel_labels.find(nidx);

        if (it != voxel_labels.end() && it->second == 0)
        {
          it->second = current_label;
          stack.push_back(nidx);
          prog->tick();
        }
      }
    }

    prog->update();
    if (prog->check_interrupt()) return {};
  }

  prog->finalize();

  std::vector<int> result(npoints);
  for (int i = 0; i < npoints; ++i)
    result[i] = voxel_labels[point_to_voxel[i]];

  return result;
}

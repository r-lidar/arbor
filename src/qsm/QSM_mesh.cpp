#include "QSM.h"
#include <algorithm>

namespace arbor::qsm {

static inline std::array<double,3> normalize(const std::array<double,3>& v)
{
  double n = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  return {v[0]/n, v[1]/n, v[2]/n};
}

static inline std::array<double,3> cross(const std::array<double,3>& a, const std::array<double,3>& b)
{
  return {
    a[1]*b[2] - a[2]*b[1],
    a[2]*b[0] - a[0]*b[2],
    a[0]*b[1] - a[1]*b[0]
  };
}

void QSM::mesh(std::vector<std::array<double,3>>& vertices, int resolution) const
{
  vertices.clear();

  // ---- sort cylinder IDs ----
  std::vector<int> ids;
  ids.reserve(cylinders_.size());
  for (const auto& kv : cylinders_)  ids.push_back(kv.first);
  std::sort(ids.begin(), ids.end());

  for (int cid : ids)
  {
    const QSMcylinder& c = cylinders_.at(cid);

    std::array<double,3> p0{c.startX, c.startY, c.startZ};
    std::array<double,3> p1{c.endX,   c.endY,   c.endZ};

    std::array<double,3> axis{p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
    double len = std::sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);

    if (len < 1e-12)  continue;

    auto dir = normalize(axis);

    // robust orthonormal basis
    std::array<double,3> z_axis{0,0,1};
    std::array<double,3> ortho1, ortho2;
    if (std::abs(dir[0]-0)<1e-6 && std::abs(dir[1]-0)<1e-6 && std::abs(dir[2]-1)<1e-6) {
      ortho1 = {1,0,0};
    } else {
      ortho1 = normalize(cross(z_axis, dir));
    }
    ortho2 = normalize(cross(dir, ortho1));

    std::vector<std::array<double,3>> bottom, top;
    bottom.reserve(resolution);
    top.reserve(resolution);

    for (int j=0; j<resolution; j++)
    {
      double theta = 2.0 * M_PI * j / resolution;
      double ct = std::cos(theta);
      double st = std::sin(theta);

      std::array<double,3> radial{
        ct * ortho1[0] + st * ortho2[0],
        ct * ortho1[1] + st * ortho2[1],
        ct * ortho1[2] + st * ortho2[2]
      };

      bottom.push_back({
        p0[0] + c.radius * radial[0],
        p0[1] + c.radius * radial[1],
        p0[2] + c.radius * radial[2]
      });

      top.push_back({
        p1[0] + c.radius * radial[0],
        p1[1] + c.radius * radial[1],
        p1[2] + c.radius * radial[2]
      });
    }

    for (auto& v : bottom) vertices.push_back(v);
    for (auto& v : top)    vertices.push_back(v);
  }
}

void QSM::tmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,3>>& faces, int resolution) const
{
  mesh(vertices, resolution);

  faces.clear();
  int num_cylinders = vertices.size() / (2 * resolution);

  for (int i = 0; i < num_cylinders; i++)
  {
    int v_offset = i * 2 * resolution;

    for (int j=0; j<resolution; j++)
    {
      int a = v_offset + j;
      int b = v_offset + ((j+1) % resolution);
      int c0 = a + resolution;
      int d = b + resolution;

      faces.push_back({a, b, d});
      faces.push_back({a, d, c0});
    }
  }
}

void QSM::qmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,4>>& faces, int resolution) const
{
  mesh(vertices, resolution);

  faces.clear();
  int num_cylinders = vertices.size() / (2 * resolution);

  for (int i = 0; i < num_cylinders; i++)
  {
    int v_offset = i * 2 * resolution;

    for (int j=0; j<resolution; j++)
    {
      int a = v_offset + j;
      int b = v_offset + ((j+1) % resolution);
      int c = b + resolution;
      int d = a + resolution;

      faces.push_back({a, b, c, d});
    }
  }
}

}

#include "QSM.h"
#include <algorithm>

namespace arbor::qsm {

void QSM::tmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,3>>& faces, MeshMode mode, int resolution) const
{
  if (mode == MeshMode::Continuous)
    tmesh_continuous(vertices, faces, resolution);

  if (mode == MeshMode::Cylinders)
    tmesh_cylinder(vertices, faces, resolution);
}

void QSM::qmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,4>>& faces, MeshMode mode, int resolution) const
{
  if (mode == MeshMode::Continuous)
    qmesh_continuous(vertices, faces, resolution);

  if (mode == MeshMode::Cylinders)
    qmesh_cylinder(vertices, faces, resolution);
}

}

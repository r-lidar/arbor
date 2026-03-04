#include "QSMbuilder.h"

#include <algorithm>

// Compute distance from point b to line ac
static double dist2line(const std::array<double, 3>& b, const std::array<double, 3>& a, const std::array<double, 3>& c)
{
  std::array<double, 3> ab, ac;
  for (int i = 0; i < 3; ++i)
  {
    ab[i] = b[i] - a[i];
    ac[i] = c[i] - a[i];
  }

  double t_num = 0.0, t_den = 0.0;
  for (int i = 0; i < 3; ++i)
  {
    t_num += ab[i] * ac[i];
    t_den += ac[i] * ac[i];
  }

  double t = t_num / t_den;
  std::array<double, 3> proj;
  for (int i = 0; i < 3; ++i)
  {
    proj[i] = a[i] + t * ac[i];
  }

  double d2 = 0.0;
  for (int i = 0; i < 3; ++i)
  {
    double diff = b[i] - proj[i];
    d2 += diff * diff;
  }

  return std::sqrt(d2);
}

void QSMbuilder::smooth_skeleton(int niter, double th)
{
  if (qsm.cylinders_.empty()) return;

  logger("Smoothing skeleton");

  // Build axis map: axis_ID -> ordered list of cyl_IDs
  std::unordered_map<int, std::vector<int>> axis_map;
  for (const auto& [id, c] : qsm)
    axis_map[c.axis_ID].push_back(id);

  // IMPORTANT: sort each axis by cyl_ID
  for (auto& [axis, vec] : axis_map) { std::sort(vec.begin(), vec.end()); }

  // Iterative smoothing
  for (int iter = 0; iter < niter; ++iter)
  {
    for (const auto& [axis_id, indices] : axis_map)
    {
      if (indices.size() < 2) continue;

      for (size_t j = 1; j < indices.size(); ++j)
      {
        int prev_id = indices[j - 1];
        int curr_id = indices[j];

        auto& prev = qsm.cylinders_[prev_id];
        auto& curr = qsm.cylinders_[curr_id];

        std::array<double,3> a = { curr.endX,  curr.endY,  curr.endZ  };
        std::array<double,3> b = { prev.startX, prev.startY, prev.startZ };
        std::array<double,3> c = { prev.endX,   prev.endY,   prev.endZ   };

        double d = dist2line(b, a, c);
        if (d <= th) continue;

        // New midpoint coordinate
        std::array<double,3> mid = {
          0.5 * (prev.startX + curr.endX),
          0.5 * (prev.startY + curr.endY),
          0.5 * (prev.startZ + curr.endZ)
        };

        // Apply displacement
        prev.endX = mid[0];
        prev.endY = mid[1];
        prev.endZ = mid[2];

        curr.startX = mid[0];
        curr.startY = mid[1];
        curr.startZ = mid[2];

        // Move all children who begin at prev.end
        auto it = qsm.children_map_.find(prev_id);
        if (it != qsm.children_map_.end())
        {
          for (int child_id : it->second)
          {
            auto& child = qsm.cylinders_[child_id];
            child.startX = mid[0];
            child.startY = mid[1];
            child.startZ = mid[2];
          }
        }
      }
    }
  }
}

#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <cmath>

using namespace Rcpp;

// Compute distance from point b to line ac
double dist2line(const std::array<double, 3>& b, const std::array<double, 3>& a, const std::array<double, 3>& c)
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

List cpp_smooth_skeleton(DataFrame qsm, int niter = 1, double th = 0)
{
  std::vector<double> startX = as<std::vector<double>>(qsm["startX"]);
  std::vector<double> startY = as<std::vector<double>>(qsm["startY"]);
  std::vector<double> startZ = as<std::vector<double>>(qsm["startZ"]);
  std::vector<double> endX   = as<std::vector<double>>(qsm["endX"]);
  std::vector<double> endY   = as<std::vector<double>>(qsm["endY"]);
  std::vector<double> endZ   = as<std::vector<double>>(qsm["endZ"]);

  std::vector<int> parent_ID = as<std::vector<int>>(qsm["parent_ID"]);
  std::vector<int> cyl_ID    = as<std::vector<int>>(qsm["cyl_ID"]);
  std::vector<int> axis_ID   = as<std::vector<int>>(qsm["axis_ID"]);

  int n = startX.size();

  std::unordered_map<int, std::vector<int>> axis_map;
  for (int i = 0; i < n; ++i)
  {
    axis_map[axis_ID[i]].push_back(i);
  }

  std::unordered_map<int, std::vector<int>> child_map;
  for (int i = 0; i < n; ++i)
  {
    child_map[parent_ID[i]].push_back(i);
  }

  for (int iter = 0; iter < niter; ++iter)
  {
    for (const auto& [axis_id, indices] : axis_map)
    {
      if (indices.size() < 2) continue;

      for (size_t j = 1; j < indices.size(); ++j)
      {
        int idx_prev = indices[j - 1];
        int idx_curr = indices[j];

        std::array<double, 3> a = {endX[idx_curr], endY[idx_curr], endZ[idx_curr]};
        std::array<double, 3> b = {startX[idx_prev], startY[idx_prev], startZ[idx_prev]};
        std::array<double, 3> c = {endX[idx_prev], endY[idx_prev], endZ[idx_prev]};

        double d = dist2line(b, a, c);

        if (d > th)
        {
          std::array<double, 3> new_coord = {
            0.5 * (startX[idx_prev] + endX[idx_curr]),
            0.5 * (startY[idx_prev] + endY[idx_curr]),
            0.5 * (startZ[idx_prev] + endZ[idx_curr])
          };

          endX[idx_prev] = new_coord[0];
          endY[idx_prev] = new_coord[1];
          endZ[idx_prev] = new_coord[2];

          startX[idx_curr] = new_coord[0];
          startY[idx_curr] = new_coord[1];
          startZ[idx_curr] = new_coord[2];

          for (int child : child_map[cyl_ID[idx_prev]])
          {
            startX[child] = new_coord[0];
            startY[child] = new_coord[1];
            startZ[child] = new_coord[2];
          }
        }
      }
    }
  }

  return List::create(
    Named("startX") = startX,
    Named("startY") = startY,
    Named("startZ") = startZ,
    Named("endX")   = endX,
    Named("endY")   = endY,
    Named("endZ")   = endZ
  );
}

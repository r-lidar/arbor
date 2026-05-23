/**
 * @file anisotopy.cpp
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

#include <vector>

#include "myomp.h"
#include "PointCloud.h"
#include "services.h"

namespace arbor::utils {

inline void eigenvalues_sym_3x3(double a, double b, double c, double d, double e, double f, double& lmin, double& lmid, double& lmax)
{
  // Matrix:
  // [ a b c ]
  // [ b d e ]
  // [ c e f ]
  const double m = (a + d + f) / 3.0;

  const double a11 = a - m;
  const double a22 = d - m;
  const double a33 = f - m;

  const double p2 = a11*a11 + a22*a22 + a33*a33 + 2.0*(b*b + c*c + e*e);

  if (p2 < 1e-30)
  {
    lmin = lmid = lmax = m;
    return;
  }

  const double p = std::sqrt(p2 / 6.0);

  // B = (1 / p) * (A - mI)
  const double b11 = a11 / p;
  const double b12 = b   / p;
  const double b13 = c   / p;
  const double b22 = a22 / p;
  const double b23 = e   / p;
  const double b33 = a33 / p;

  const double detB = b11*(b22*b33 - b23*b23) - b12*(b12*b33 - b23*b13) + b13*(b12*b23 - b22*b13);

  double r = detB * 0.5;
  r = std::max(-1.0, std::min(1.0, r));

  const double phi = std::acos(r) / 3.0;
  const double two_p = 2.0 * p;

  lmax = m + two_p * std::cos(phi);
  lmin = m + two_p * std::cos(phi + 2.0*M_PI/3.0);
  lmid = 3.0*m - lmax - lmin;
}

std::vector<float> anisotropy(const PointCloud& cloud, int k)
{
  const int n = cloud.size();

  cloud.build_index();

  std::vector<float> out(n);


  auto pb = ServiceLocator::make_progress(n, "Anisotropy");
  std::atomic<bool> abort(false);

  #pragma omp parallel
  {
    std::vector<unsigned int> idx(k);
    std::vector<double> dist(k);
    double q[3];

    #pragma omp for schedule(static)
    for (int i = 0; i < n; ++i)
    {
      if (abort.load(std::memory_order_relaxed)) continue;
      if(pb->check_interrupt()) abort = true;
      pb->tick();

      cloud.get_point(i, q);
      cloud.knn(q, k, idx, dist);

      // ---- mean ----
      double mx=0, my=0, mz=0;
      for (int j=0; j<k; ++j)
      {
        mx += cloud.get_x(idx[j]);
        my += cloud.get_y(idx[j]);
        mz += cloud.get_z(idx[j]);
      }
      mx /= k; my /= k; mz /= k;

      // ---- covariance ----
      double cxx=0, cxy=0, cxz=0, cyy=0, cyz=0, czz=0;
      for (int j=0; j<k; ++j)
      {
        double dx = cloud.get_x(idx[j]) - mx;
        double dy = cloud.get_y(idx[j]) - my;
        double dz = cloud.get_z(idx[j]) - mz;

        cxx += dx*dx;
        cxy += dx*dy;
        cxz += dx*dz;
        cyy += dy*dy;
        cyz += dy*dz;
        czz += dz*dz;
      }

      const double inv = 1.0 / (k - 1);
      cxx *= inv; cxy*=inv; cxz*=inv;
      cyy *= inv; cyz*=inv; czz*=inv;

      const double trace = cxx + cyy + czz;
      if (trace < 1e-14)
      {
        out[i] = 0.0;
        continue;
      }

      double lmin, lmid, lmax;
      eigenvalues_sym_3x3(cxx, cxy, cxz, cyy, cyz, czz, lmin, lmid, lmax);

      out[i] = static_cast<float>((lmax > 1e-14) ? (lmax - lmin) / lmax : 0.0);
    }

    pb->tick();
  }

  if (abort.load()) throw std::runtime_error("Computation aborted");

  return out;
}

}

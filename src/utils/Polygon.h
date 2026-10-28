/**
 * @file Polygon.h
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

#ifndef POLYGON_H
#define POLYGON_H

#include <vector>
#include <cmath>
#include <cstddef>
#include <limits>

namespace arbor::utils
{

// A simple 2D polygon defined by a single ring (no holes, no multipolygon).
// Points lying on the boundary are considered inside (same as sf::st_intersects).
class Polygon
{
public:
  struct Point { double x, y; };
  using Ring = std::vector<Point>;

  Polygon() = default;

  // The ring can be open or closed (first == last). A ring with less than 3 distinct vertices
  // gives an empty polygon.
  explicit Polygon(const Ring& ring);

  bool empty() const { return vertices.empty(); }
  const Ring& ring() const { return vertices; }

  bool point_in_polygon(double x, double y) const;

  // Shrunk polygon (negative buffer). Only convex polygons are supported. The result contains the
  // points of the polygon that are at least d >= 0 from its boundary. It is empty if d is too large
  // for the polygon (same as sf::st_buffer() with distance -d).
  Polygon shrink(double d) const;

  bool is_convex() const;

  // Convex hull (Andrew's monotone chain). Returns an empty polygon if the hull is degenerate.
  static Polygon convex_hull(std::vector<Point> pts);

  // Parallel convex hull of n points accessed via get_x(i) and get_y(i). Points that are strictly
  // inside the octagon defined by the 8 extreme points (Akl-Toussaint heuristic) are discarded
  // in parallel and the hull is computed on the remaining points. Non finite points are ignored.
  template <typename GetX, typename GetY>
  static Polygon convex_hull(size_t n, GetX get_x, GetY get_y);

private:
  Ring vertices;
  double xmin = 0, ymin = 0, xmax = 0, ymax = 0;
};

template <typename GetX, typename GetY>
Polygon Polygon::convex_hull(size_t n, GetX get_x, GetY get_y)
{
  constexpr double inf = std::numeric_limits<double>::infinity();

  // 1. Find the 8 extreme points: min y, max x-y, max x, max x+y, max y, min x-y, min x, min x+y
  // These points are on the hull and are listed in counter-clockwise order.
  struct Extremes
  {
    double val[8] = {inf, -inf, -inf, -inf, -inf, inf, inf, inf};
    Point  pt[8];
  };

  auto update = [](Extremes& e, double x, double y)
  {
    const double v[8] = {y, x - y, x, x + y, y, x - y, x, x + y};
    for (int k = 0; k < 8; ++k)
    {
      bool is_min = (k == 0 || k >= 5);
      if ((is_min && v[k] < e.val[k]) || (!is_min && v[k] > e.val[k]))
      {
        e.val[k] = v[k];
        e.pt[k] = {x, y};
      }
    }
  };

  Extremes global;

  #pragma omp parallel
  {
    Extremes local;

    #pragma omp for schedule(static)
    for (size_t i = 0; i < n; ++i)
    {
      double x = get_x(i);
      double y = get_y(i);
      if (!std::isfinite(x) || !std::isfinite(y)) continue;
      update(local, x, y);
    }

    #pragma omp critical
    {
      for (int k = 0; k < 8; ++k)
      {
        if (std::isfinite(local.val[k])) update(global, local.pt[k].x, local.pt[k].y);
      }
    }
  }

  if (!std::isfinite(global.val[0])) return Polygon(); // no finite point

  // 2. Octagon made of the distinct extreme points (counter-clockwise)
  Ring oct;
  for (int k = 0; k < 8; ++k)
  {
    const Point& p = global.pt[k];
    if (oct.empty() || oct.back().x != p.x || oct.back().y != p.y) oct.push_back(p);
  }
  while (oct.size() > 1 && oct.front().x == oct.back().x && oct.front().y == oct.back().y) oct.pop_back();
  bool filter = oct.size() >= 3;

  // 3. Keep only points that are not strictly inside the octagon
  std::vector<Point> candidates;

  #pragma omp parallel
  {
    std::vector<Point> local;

    #pragma omp for schedule(static) nowait
    for (size_t i = 0; i < n; ++i)
    {
      double x = get_x(i);
      double y = get_y(i);
      if (!std::isfinite(x) || !std::isfinite(y)) continue;

      bool strictly_inside = filter;
      for (size_t k = 0; filter && k < oct.size(); ++k)
      {
        const Point& a = oct[k];
        const Point& b = oct[(k + 1) % oct.size()];
        if ((b.x - a.x) * (y - a.y) - (b.y - a.y) * (x - a.x) <= 0)
        {
          strictly_inside = false;
          break;
        }
      }

      if (!strictly_inside) local.push_back({x, y});
    }

    #pragma omp critical
    candidates.insert(candidates.end(), local.begin(), local.end());
  }

  return convex_hull(std::move(candidates));
}

}

#endif

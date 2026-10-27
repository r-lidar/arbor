/**
 * @file Polygon.cpp
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

#include "Polygon.h"

#include <algorithm>

Polygon::Polygon(const Ring& ring)
{
  vertices.reserve(ring.size());
  for (const Point& p : ring)
  {
    if (!std::isfinite(p.x) || !std::isfinite(p.y)) continue;
    if (vertices.empty() || vertices.back().x != p.x || vertices.back().y != p.y) vertices.push_back(p);
  }

  // Open the ring if closed
  while (vertices.size() > 1 && vertices.front().x == vertices.back().x && vertices.front().y == vertices.back().y) vertices.pop_back();

  if (vertices.size() < 3)
  {
    vertices.clear();
    return;
  }

  xmin = xmax = vertices[0].x;
  ymin = ymax = vertices[0].y;
  for (const Point& p : vertices)
  {
    xmin = std::min(xmin, p.x);
    xmax = std::max(xmax, p.x);
    ymin = std::min(ymin, p.y);
    ymax = std::max(ymax, p.y);
  }
}

// Even-odd ray casting. Points on an edge are inside.
bool Polygon::point_in_polygon(double x, double y) const
{
  if (empty()) return false;
  if (x < xmin || x > xmax || y < ymin || y > ymax) return false;

  bool inside = false;
  size_t n = vertices.size();
  for (size_t i = 0, j = n - 1; i < n; j = i++)
  {
    const Point& a = vertices[i];
    const Point& b = vertices[j];

    // Exactly on the segment [a, b]
    double cross = (b.x - a.x) * (y - a.y) - (b.y - a.y) * (x - a.x);
    if (cross == 0 &&
        x >= std::min(a.x, b.x) && x <= std::max(a.x, b.x) &&
        y >= std::min(a.y, b.y) && y <= std::max(a.y, b.y))
      return true;

    if ((a.y > y) != (b.y > y))
    {
      double xi = a.x + (y - a.y) * (b.x - a.x) / (b.y - a.y);
      if (x < xi) inside = !inside;
    }
  }
  return inside;
}

double Polygon::distance_to_boundary(double x, double y) const
{
  double dmin = std::numeric_limits<double>::infinity();
  size_t n = vertices.size();
  for (size_t i = 0, j = n - 1; i < n; j = i++)
  {
    const Point& a = vertices[j];
    const Point& b = vertices[i];
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double len2 = dx * dx + dy * dy;
    double t = (len2 > 0) ? ((x - a.x) * dx + (y - a.y) * dy) / len2 : 0;
    t = std::clamp(t, 0.0, 1.0);
    double px = a.x + t * dx - x;
    double py = a.y + t * dy - y;
    dmin = std::min(dmin, px * px + py * py);
  }
  return std::sqrt(dmin);
}

bool Polygon::contains(double x, double y, double margin) const
{
  if (!point_in_polygon(x, y)) return false;
  if (margin <= 0) return true;
  return distance_to_boundary(x, y) >= margin;
}

Polygon Polygon::convex_hull(std::vector<Point> pts)
{
  std::sort(pts.begin(), pts.end(), [](const Point& a, const Point& b) { return a.x < b.x || (a.x == b.x && a.y < b.y); });
  pts.erase(std::unique(pts.begin(), pts.end(), [](const Point& a, const Point& b) { return a.x == b.x && a.y == b.y; }), pts.end());

  if (pts.size() < 3) return Polygon();

  auto cross = [](const Point& o, const Point& a, const Point& b) { return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x); };

  // Andrew's monotone chain. Collinear points are removed. Result is counter-clockwise.
  Ring hull(2 * pts.size());
  size_t k = 0;
  for (size_t i = 0; i < pts.size(); ++i)
  {
    while (k >= 2 && cross(hull[k - 2], hull[k - 1], pts[i]) <= 0) k--;
    hull[k++] = pts[i];
  }
  for (size_t i = pts.size() - 1, t = k + 1; i > 0; --i)
  {
    while (k >= t && cross(hull[k - 2], hull[k - 1], pts[i - 1]) <= 0) k--;
    hull[k++] = pts[i - 1];
  }
  hull.resize(k - 1);

  return Polygon(hull);
}

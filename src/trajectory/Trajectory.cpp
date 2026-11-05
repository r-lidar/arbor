/**
 * @file Trajectory.cpp
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

#include "Trajectory.h"

#include <algorithm>

namespace arbor {

Trajectory::Trajectory(std::vector<Pose> poses)
{
  // Files are written in time order: the sort is then a check. Stable, so that positions at the same time keep the
  // order of the file.
  auto by_time = [](const Pose &a, const Pose &b) { return a.time < b.time; };
  if (!std::is_sorted(poses.begin(), poses.end(), by_time))
    std::stable_sort(poses.begin(), poses.end(), by_time);

  times.reserve(poses.size());
  positions.reserve(poses.size());
  for (const Pose &pose : poses)
  {
    times.push_back(pose.time);
    positions.push_back(pose.position);
  }
}

// ── Queries ─────────────────────────────────────────────────────────────────────

Trajectory::Point Trajectory::position_at(double t) const
{
  if (times.size() < 2)
    return empty() ? Point{0, 0, 0} : positions.front(); // no segment

  return interpolate(segment_at(t), t);
}

// The previous segment, or the next one (t moved past its end), before a binary search: queries in time order stay
// in O(1).
Trajectory::Point Trajectory::position_at(double t, size_t &hint) const
{
  if (times.size() < 2)
    return empty() ? Point{0, 0, 0} : positions.front(); // no segment

  const size_t segment_count = times.size() - 1;
  auto is_in = [this, t](size_t segment) { return times[segment] <= t && t < times[segment + 1]; };

  if (hint < segment_count && is_in(hint))
    return interpolate(hint, t);

  if (hint + 1 < segment_count && is_in(hint + 1))
  {
    ++hint;
    return interpolate(hint, t);
  }

  hint = segment_at(t);
  return interpolate(hint, t);
}

// ── Change ──────────────────────────────────────────────────────────────────────

void Trajectory::translate(double dx, double dy, double dz)
{
  for (Point &p : positions)
  {
    p.x += dx;
    p.y += dy;
    p.z += dz;
  }
}

// ── Helpers ─────────────────────────────────────────────────────────────────────

// The last position at or before t starts the segment (upper_bound: after the positions at the same time, so the last
// of them is used). Before the first position: segment 0; at or after the last one: the last segment.
size_t Trajectory::segment_at(double t) const
{
  const size_t index = static_cast<size_t>(std::upper_bound(times.begin(), times.end(), t) - times.begin());
  const size_t last  = times.size() - 2;
  return index == 0 ? 0 : std::min(index - 1, last);
}

// The fraction of the segment at t, clamped to [0, 1] (queries outside the recorded times get an end). A segment of
// null duration (two positions at the same time, the last segment): its end, the last position at that time.
Trajectory::Point Trajectory::interpolate(size_t segment, double t) const
{
  const double start    = times[segment];
  const double duration = times[segment + 1] - start;
  if (duration <= 0.0)
    return positions[segment + 1];

  const double fraction = std::clamp((t - start) / duration, 0.0, 1.0);
  const Point &from     = positions[segment];
  const Point &to       = positions[segment + 1];
  return {from.x + fraction * (to.x - from.x),
          from.y + fraction * (to.y - from.y),
          from.z + fraction * (to.z - from.z)};
}

} // namespace arbor

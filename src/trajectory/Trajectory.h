/**
 * @file Trajectory.h
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

#ifndef ARBOR_TRAJECTORY_H
#define ARBOR_TRAJECTORY_H

#include <cstddef>
#include <vector>

namespace arbor {

// Positions of the scanner ordered by time, to answer "where was it at time t" fast, in a hot loop: the position
// between the two recorded positions around t, linearly interpolated.
//
// Stored as two parallel vectors sorted by time (times, positions): the binary search of a query walks the times
// only, contiguous doubles (8 B per position), and never allocates. position_at(t) costs O(log n);
// position_at(t, hint), for queries that follow each other in time, O(1) while t stays in or next to the segment
// of the previous query.
//
// Times are in the unit of the file (GPS seconds...), positions in any frame (double precision: world coordinates
// fit). Several positions at the same time are kept; a query at that time gets the last of them.
class Trajectory
{
public:
  struct Point { double x, y, z; };

  // One recorded position.
  struct Pose
  {
    double time;
    Point  position;
  };

  Trajectory() = default;
  explicit Trajectory(std::vector<Pose> poses); // any order: sorted by time (stable: equal times keep their order)

  // Recorded positions, by increasing time
  bool                       empty()           const { return times.empty(); }
  size_t                     size()            const { return times.size(); }
  const std::vector<double> &get_times()       const { return times; }
  const std::vector<Point>  &get_positions()   const { return positions; }
  double                     get_start_time()  const { return times.front(); } // not empty() only
  double                     get_end_time()    const { return times.back(); }  // not empty() only
  bool                       covers(double t)  const { return !empty() && t >= times.front() && t <= times.back(); }

  // Position at time t, linearly interpolated between the recorded positions around it. Clamped outside the recorded
  // times (the first or the last position: check covers() when it matters). Empty: (0, 0, 0).
  Point position_at(double t) const;

  // Same, for a sequence of queries: `hint` is the segment of the previous query (start at 0, keep the variable
  // between the calls of one thread). The segment i joins the positions i and i + 1.
  Point position_at(double t, size_t &hint) const;

  // Moves every position (e.g. from world coordinates to a local frame).
  void translate(double dx, double dy, double dz);

private:
  size_t segment_at(double t) const;                   // binary search: the segment i with times[i] <= t < times[i + 1]
  Point  interpolate(size_t segment, double t) const;  // along the segment, clamped to its ends

  std::vector<double> times;     // increasing
  std::vector<Point>  positions; // positions[i] at times[i]
};

} // namespace arbor

#endif // ARBOR_TRAJECTORY_H

/**
 * @file fitting_strategy.h
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

#pragma once

#include "fitting_types.h"
#include <cstdint>
#include <vector>

namespace arbor::utils::fitting {

// ---------------------------------------------------------------------------
// Abstract base for all 2-D shape-fitting strategies.
// Concrete implementations receive a cloud of Vec3 points that have already
// been projected onto the XY plane (Z is ignored by the fitting math but
// preserved for the output center / contour).
// ---------------------------------------------------------------------------
class ShapeFitter
{
public:
  virtual ~ShapeFitter() = default;
  virtual FittingResult fit(const std::vector<Vec3>& points, double tolerance) = 0;

  /**
   * @brief Measure how many angular degrees around the fitted center
   *        are covered by inlier points.
   *
   * The circle is divided into 10° bins; the method counts how many
   * consecutive bins are occupied.  Used by every fitting strategy to
   * populate FittingResult::arc_coverage_deg.
   */
  static double calculate_arc_coverage(const std::vector<Vec3>& points,
                                        const std::vector<int>&  inliers,
                                        const Vec3&              center);

protected:
  double m_zmean = 0.0;   // mean Z of the input cloud, set at the top of each fit()
};

// ---------------------------------------------------------------------------
// Bitmask of fitting algorithms to enable.
// Combine with | and pass to CrossSectionFitter::fit().
// ---------------------------------------------------------------------------
enum class FitMode : uint32_t
{
  None         = 0,
  Circle       = 1 << 0,   // 0x01 — single circle (RANSAC)
  Ellipse      = 1 << 1,   // 0x02 — ellipse (algebraic + RANSAC)
  Fourier5     = 1 << 2,   // 0x04 — Fourier contour, 5 harmonics (complex stem)
  MultiCircle2 = 1 << 3,   // 0x08 — two circles (double stem)
  MultiCircle3 = 1 << 4,   // 0x10 — three circles (triple stem)
  Fourier10    = 1 << 5,   // 0x20 — Fourier contour, 10 harmonics (tropical buttress)

  // Convenience presets
  Basic    = Circle | Ellipse,
  Standard = Basic  | Fourier5,
  Buttress = Standard | Fourier10,
  Full     = Standard | MultiCircle2 | MultiCircle3 | Fourier10,
};

inline FitMode operator|(FitMode a, FitMode b)
{
  return static_cast<FitMode>(static_cast<uint32_t>(a) | static_cast<uint32_t>(b));
}

inline FitMode operator&(FitMode a, FitMode b)
{
  return static_cast<FitMode>(static_cast<uint32_t>(a) & static_cast<uint32_t>(b));
}

inline bool has(FitMode flags, FitMode flag)
{
  return (flags & flag) != FitMode::None;
}

} // namespace arbor::utils::fitting

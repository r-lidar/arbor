/**
 * @file FitQuality.h
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

#ifndef FIT_QUALITY_H
#define FIT_QUALITY_H

#include "fitting.h"
#include <algorithm>

namespace arbor::utils::fitting {

/**
 * @brief Named acceptance criteria for a single circle-fitting context.
 *
 *
 * Usage:
 *   auto res = fitter.fit(quality.ransac_thr);
 *   if (quality.accept(res)) { ... }                    // no prior radius
 *   if (quality.accept(res, einfo.data.radius)) { ... } // with shrinkage guard
 */
struct FitQuality
{
  float min_radius  = 0.0f;   // Minimum acceptable fitted radius (m).
  float min_arc_deg = 0.0f;   // Minimum arc coverage required (degrees).
  float min_inlier  = 0.0f;   // Minimum inlier percentage required (%).
  float max_shrink  = 1.0f;   // Maximum fractional shrinkage allowed relative
  float ransac_tolerance = 0.04f;  // RANSAC inlier-distance threshold (m) passed

  // -----------------------------------------------------------------------
  // Acceptance predicates
  // -----------------------------------------------------------------------

  /**
   * Accept a fit result with no prior-radius context.
   * Tests radius, arc coverage, and inlier percentage only.
   */
  bool accept(const utils::fitting::FittingResult& res) const noexcept
  {
    return
    res.success             == true
    && res.radius           >= min_radius
    && res.arc_coverage_deg  > min_arc_deg
    && res.inlier_percentage > min_inlier;
  }

  /**
   * Accept a fit result with a shrinkage guard against a known prior radius.
   *
   * The fit is rejected when the new radius has shrunk by more than
   * `max_shrink` relative to `current_radius`:
   *   ratio = (res.radius - current_radius) / current_radius
   *   ratio > -max_shrink  must hold
   *
   * The guard is skipped when max_shrink == 1.0f (default) or when
   * current_radius <= 0.
   */
  bool accept(const utils::fitting::FittingResult& res, float current_radius) const noexcept
  {
    if (!accept(res)) return false;
    if (max_shrink < 1.0f && current_radius > 0.0f)
    {
      const float ratio = (res.radius - current_radius) / current_radius;
      if (ratio < -max_shrink) return false;
    }
    return true;
  }

  static float quality_score(const utils::fitting::FittingResult& res) noexcept
  {
    const float arc_fraction    = std::min(res.arc_coverage_deg, 360.0f) / 360.0f;
    const float inlier_fraction = std::min(res.inlier_percentage, 100.0f) / 100.0f;
    return arc_fraction * inlier_fraction;
  }

  // -----------------------------------------------------------------------
  // Named presets  (one per calling context)
  // -----------------------------------------------------------------------

  /**
   * First measurement pass - measure_radii().
   *
   * Permissive arc threshold (180°) so partial scans are still used.
   * Minimum radius (0.03 m) rejects fits on very thin twigs whose
   * conic-allometry pre-filter already left them unset.
   */
  static constexpr FitQuality standard_preset() noexcept
  {
    return { /*min_radius*/   0.03f,
             /*min_arc_deg*/ 180.0f,
             /*min_inlier*/   30.0f,
             /*max_shrink*/    1.0f};
  }

  static constexpr FitQuality low_preset() noexcept
  {
    return { /*min_radius*/   0.00f,
             /*min_arc_deg*/ 120.0f,
             /*min_inlier*/   50.0f,
             /*max_shrink*/    1.0f};
  }

  /**
   * Second-pass refinement - refine_radii() and refine_radii_broken().
   *
   * Requires near-complete arc coverage (300°) and high inlier fraction
   * (70%) so only well-observed cylinders are updated.  The shrinkage cap
   * of 10% prevents a noisy fit from drastically reducing an already-good
   * radius.
   *
   * Note: refine_radii() adds an independent pre-condition
   *   `einfo.data.radius >= 0.04`
   * that guards whether fitting is attempted at all; that threshold is
   * intentionally kept at the call site because it concerns the state of
   * the cylinder, not the quality of the fit result.
   */
  static constexpr FitQuality accurate_preset() noexcept
  {
    return { /*min_radius*/    0.0f,
             /*min_arc_deg*/ 300.0f,
             /*min_inlier*/   70.0f,
             /*max_shrink*/    0.1f};
  }

  /**
   * Broken-tree ring detection - build_skeleton().
   *
   * High arc and inlier thresholds (same as refinement) to ensure only
   * clean, well-sampled cross-sections increment the broken-tree counter.
   * Minimum radius (0.06 m) excludes small-branch clusters that happen to
   * be round from triggering the broken-tree heuristic.
   */
  static constexpr FitQuality accurate_and_large_preset() noexcept
  {
    return { /*min_radius*/   0.06f,
             /*min_arc_deg*/ 300.0f,
             /*min_inlier*/   70.0f,
             /*max_shrink*/    1.0f};
  }
};

} // namespace arbor::qsm

#endif // FIT_QUALITY_H

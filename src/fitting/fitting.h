/**
 * @file fitting.h
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

/**
 * Convenience umbrella header — includes every fitting component.
 * Existing callers that include only "fitting.h" continue to work unchanged.
 * For finer-grained dependencies, include the individual headers directly.
 */

#pragma once

#include "fitting_types.h"
#include "fitting_circle_math.h"
#include "fitting_strategy.h"
#include "fitting_circle.h"
#include "fitting_ellipse.h"
#include "fitting_fourier.h"
#include "fitting_multicircle.h"
#include "fitting_orbicular.h"

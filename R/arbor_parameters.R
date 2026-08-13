# @file arbor_parameters.R
# Project: Arbor
#
# Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' Parameters
#'
#' The built-in internal variable `arbor_parameters_default` contains all the parameters used
#' by the arbor functions. The variable itself is a nested list that is intentionally not documented.
#' There are numerous parameters controlling thresholds, neighborhoods, clustering, noise filtering,
#' and pathfinder configuration. If users had to adjust and understand all of them, the software would
#' quickly become unusable. While every parameter can be configured, we deliberately avoid documenting
#' them to prevent endless tweaking. The pipeline MUST work reliably out of the box with the default
#' parameters. See the [Arbor book](<placeholder>)
#' @rdname parameters
#' @name parameters
#' @examples
#' params = arbor_parameters_default
#' params$global$cut_above_ground <- 0.25
#' @export
arbor_parameters_default = list()

#' Constants
#'
#' Built-in variables. See the [Arbor book](<placeholder>) for more details.
#'
#' - `ARBORTREE` is used in the `UserData` attributes. It flags points that are valid trees.
#'   Usually, all points of interest are flagged as `ARBORTREE`.
#' - `ARBORLOW` is used in the `UserData` attributes. It flags low points that are not processed by
#'   Arbor at all.
#' - `ARBORUNDERSTORY` is used in the `UserData` attributes. It flags small trees that are not
#'   processed to build QSMs. See \link{flag_small_trees}.
#' - `ARBORBUFFER` is used in the `UserData` attributes. It flags points in the buffer area. See
#'   \link{flag_buffer}.
#' @export
#' @md
#' @rdname constants
ARBORTREE <- 0L

#' @export
#' @rdname constants
ARBORLOW <- 1L
#' @export
#' @rdname constants
ARBORUNDERSTORY <- 2L
#' @export
#' @rdname constants
ARBORBUFFER <- 3L

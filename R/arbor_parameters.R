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
#' Default is suitable for most dataset. Occlusion is suitable for lower quality datasets
#' @rdname parameters
#' @name parameters
#' @export
arbor_parameters_default = list()

#' @export
#' @rdname parameters
#' @name parameters
arbor_parameters_occlusion = list()

#' @export
#' @rdname parameters
ARBORTREE <- 0L
#' @export
#' @rdname parameters
ARBORLOW <- 1L
#' @export
#' @rdname parameters
ARBORUNDERSTORY <- 2L
#' @export
#' @rdname parameters
ARBORBUFFER <- 3L

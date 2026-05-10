# @file crs.R
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

#' Get or set the projection of a QSM and QSF objects
#'
#' Get or set the projection of a QSM and QSF objects
#'
#' @param x a QSM or QSF
#' @param value see sf::st_crs
#' @param ... Unused.
#'
#' @export
#' @importFrom sf st_crs
#' @importFrom sf st_crs<-
#' @name st_crs
#' @md
NULL

#' @export
#' @rdname st_crs
st_crs.qsm = function(x, ...) { attr(x, "crs") }

#' @export
#' @rdname st_crs
st_crs.qsf = function(x, ...) { attr(x, "crs") }


#' @export
#' @rdname st_crs
`st_crs<-.qsm` = function(x, value) { attr(x, "crs") = value ; x }

#' @export
#' @rdname st_crs
`st_crs<-.qsf` = function(x, value) { attr(x, "crs") = value ; x }

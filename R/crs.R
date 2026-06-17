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
st_crs.qsm = function(x, ...)
{
  crs <- attr(x, "crs")
  if (is.null(crs)) crs <- sf::NA_crs_
  else if (crs == "") crs <- sf::NA_crs_
  else crs <- sf::st_crs(crs)
}

#' @export
#' @rdname st_crs
st_crs.qsf = function(x, ...)
{
  sf::st_crs(x[[1]])
}


#' @export
#' @rdname st_crs
`st_crs<-.qsm` = function(x, value) { attr(x, "crs") = clean_crs(sf::st_crs(value)) ; x }

#' @export
#' @rdname st_crs
`st_crs<-.qsf` = function(x, value)
{
  wkt = clean_crs(sf::st_crs(value))
  lapply(x, function(y) { attr(y, "crs") = wkt ; y })
}

clean_crs <- function(crs)
{
  if (inherits(crs, "crs")) crs <- crs$wkt
  if (is.na(crs)) return("")
  crs <- gsub("\\s+", " ", crs)
  crs <- gsub("\\s*([\\[\\],])\\s*", "\\1", crs)
  trimws(crs)
}

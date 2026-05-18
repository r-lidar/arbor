# @file qsf_print.R
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

#' Print a Quantitative Structure Forest (QSF)
#'
#' Displays a short summary of a Quantitative Structure Forest followed by
#' the standard \code{data.table} print output.
#'
#' @param x A \code{qsm} object.
#' @param ... Unused (for S3 compatibility).
#' @method print qsf
#' @return Invisibly returns \code{x}.
#' @export
print.qsf <- function(x, ...)
{
  # number of cylinders
  n_tree <- length(x)
  n_cyl  <- sum(sapply(x, nrow))
  crs = attr(x, "crs")
  dbh = qsm_dbh(x)
  ba  = sum(pi*(dbh$dbh/2)^2)
  if (is.null(crs)) crs = sf::NA_crs_

  # total volume
  total_volume <- sum(sapply(x, qsm_volume))

  cat("QSF\n")
  cat(sprintf("Trees       : %d\n", n_tree))
  cat(sprintf("Cylinders   : %d\n", n_cyl))
  cat(sprintf("Basal area  : %.2f m\u00b2\n", ba))
  cat(sprintf("Volume      : %.2f m\u00b3\n", total_volume))
  cat("Coord. ref. :", crs$Name, "\n")

  invisible(x)
}

# @file qsm_print.R
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

#' Print a Quantitative Structure Model (QSM)
#'
#' Displays a short summary of a Quantitative Structure Model followed by
#' the standard \code{data.table} print output.
#'
#' @param x A \code{qsm} object.
#' @param ... Unused (for S3 compatibility).
#' @method print qsm
#' @return Invisibly returns \code{x}.
#' @export
print.qsm <- function(x, ...)
{
  n_cyl <- nrow(x)
  volume <- 0
  height <- 0
  dbh <- 0
  if (n_cyl > 0)
  {
    volume <- qsm_volume(x)
    zmin <- min(c(x$startZ, x$endZ), na.rm = TRUE)
    zmax <- max(c(x$startZ, x$endZ), na.rm = TRUE)
    height <- zmax - zmin
    dbh = qsm_dbh(x)$dbh
  }
  #crs = attr(x, "crs")
  #if (is.null(crs)) crs = sf::NA_crs_

  msg <- qsm_message(x)

  cat(sprintf("Class       : QSM\n"))
  cat(sprintf("Cylinders   : %d\n", n_cyl))
  cat(sprintf("Diameter    : %.1f cm\n", dbh*100))
  cat(sprintf("Height      : %.1f m\n", height))
  cat(sprintf("Volume      : %.2f m\u00b3\n", volume))
  #cat("Coord. ref. :", crs$Name, "\n")
  if (length(msg) > 0) cat("Message     :", msg, "\n")

  invisible(x)
}

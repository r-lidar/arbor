# @file qsm.R
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

#' Generate a QSM from a single tree point cloud
#'
#' This function processes a tree point cloud to generate a Quantitative Structure Model (QSM)
#' See the [Arbor book](https://r-lidar.github.io/arbor_book/) for mode details.
#'
#' @param tree A `LAS` object from lidR containing a single tree point cloud. Only
#' the point labelled as wood will be used for QSM.
#' @param params list See \link{parameters}.
#' @examples
#' f <- system.file("extdata", "tree_qsm.laz", package="arbor")
#' tree <- lidR::readLAS(f)
#' qsm = qsm(tree)
#' \dontrun{
#' x = plot_semantic(tree)
#' plot_qsm(qsm, add = x, color = "branch_order", cylinder = TRUE)
#' }
#' @export
#' @md
#' @seealso \link{qsm_write} \link{qsm_read} \link{qsm_dbh} \link{qsm_stats}
qsm =  function(tree, params = arbor_parameters_default)
{
  qsm <- qsm_cpp(tree@data, params)
  qsm <- qsm_finalize(qsm)
  st_crs(qsm) <- st_crs(tree)
  qsm
}

qsm_finalize = function(qsm)
{
  data.table::setDT(qsm)
  order <- c("startX", "startY", "startZ", "endX", "endY", "endZ", "cyl_ID", "parent_ID", "axis_ID", "branch_order", "dist_to_root", "subtree_length", "radius")
  data.table::setcolorder(qsm, order)
  qsm <- as_qsm(qsm)
  msg = qsm_message(qsm)
  if (length(msg) > 0) warning(msg, call. = FALSE)
  qsm
}

as_qsm <- function(x)
{
  data.table::setDT(x)
  class(x) <- c("qsm", class(x))
  x
}


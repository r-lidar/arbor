# @file segment_instance.R
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

#' Individual Tree Instance Segmentation
#'
#' This function performs individual tree instance segmentation by constructing a point-based network
#' and computing the least-cost path for each point to predefined seed points. The goal is to assign
#' each point in the point cloud a `treeID` corresponding to its nearest seed, based on the least-cost
#' path computed through the network. The process incorporates both spatial and heuristic constraints
#' to segment trees effectively. See the [Arbor book](https://r-lidar.github.io/arbor_book/) for mode details.
#'
#' @param las A `LAS` object (from the `lidR` package) containing the point cloud data. It must have
#'   an attribute `foliage`, which can be generated using prior semantic segmentation.
#' @param seeds An `sf` object containing points that represent tree seeds. Each point must have
#'   xyz coordinates and a `treeID`. Seeds can be generated using \link{find_seeds}.
#' @param params list See \link{parameters}.
#'
#' @return A `LAS` object with an additional attribute `treeID`, indicating the ID of the tree to
#'   which each point belongs. `treeID` is never `NA`: it is either a valid ID `> 0`, or `0` for
#'   points that do not belong to any tree (points below `cut_above_ground` and points that cannot
#'   be connected to any seed).
#'
#' @md
#' @export
#' @rdname segment_instance
#' @name segment_instance
#' @importFrom Rcpp sourceCpp
#' @seealso \link{find_seeds}, \link{segment_semantic}
segment_instance = function(las, seeds, params = arbor_parameters_default)
{
  stop_if_not_tls(las)
  las@data$treeID = 0L
  segment_instance_cpp(las@data, seeds@data, params)
  las <- lidR::add_lasattribute_manual(las, name = "treeID", desc = "Unique ID per tree", type = "int")
  return(las)
}

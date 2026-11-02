# @file segment_semantic.R
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

#' Wood-Foliage Semantic Segmentation
#'
#' This function performs wood-foliage semantic segmentation by constructing a network of points
#' using k-nearest neighbors (KNN). A pathfinding algorithm navigates this network to identify the
#' tree skeleton by finding the least-cost path from every point in space to the ground. Each point
#' is classified as wood or foliage based on two main criteria: (1) Proximity to the detected skeleton and
#' (2) wood likelihood exceeding a given threshold. The function \link{wood_likelihood} must be applied
#' first. The point cloud must have an attribute 'hag' (see \link[lidR:height_above_ground]{height_above_ground})
#' Finally, a connected component step removes small clusters incorrectly classified as wood.\cr\cr
#' See the [Arbor book](https://r-lidar.github.io/arbor_book/) for mode details.
#'
#' @param las A LAS object.
#' @param params list See \link{parameters}.
#'
#' @return A `LAS` object with an additional attribute `foliage`. 0 - wood, 1 and 2 foliage. The `UserData`
#' attribute is also assigned two constant: `ARBORLOW` for points below the cut threshold ans `ARBORTREE`
#' for all others points
#' @md
#' @rdname segment_semantic
#' @name segment_semantic
#' @export
segment_semantic = function(las, params = arbor_parameters_default)
{
  # Before to enter the semantic segmentation we must ensure memory is allocated
  # by R for attributes 'passage', 'foliage' and 'UserData'
  las@data$UserData <- ARBORTREE
  las@data$passage  <- 0L
  las@data$foliage  <- 1L
  las <- lidR::add_lasattribute_manual(las, name = "passage", desc = "passage points", type = "int")
  las <- lidR::add_lasattribute_manual(las, name = "foliage", desc = "foliage: 1 or 2 wood: 0", type = "char")

  segment_semantic_cpp(las@data, params)
  return(las)
}

# @file segment_post_prod.R
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

#' Flag Trees for Exclusion
#'
#' Utilities to flag trees that should be excluded from downstream processing.
#' These functions assign values to the `UserData` attribute of the point cloud in order to later
#' filter data on the fly. See the [Arbor book](https://r-lidar.github.io/arbor_book/)
#'
#' `flag_small_trees()` flags trees below a given height threshold, which are often
#' poorly segmented or correspond to understory vegetation.\cr\cr
#' `flag_buffer()` flags trees located near the edges of a point cloud by clipping
#' trees whose seed points fall outside a buffered polygon.\cr\cr
#' See the [Arbor book](https://r-lidar.github.io/arbor_book/) for more details.
#'
#' @param las A LAS object from lidR containing segmented trees.
#' @param max_height Numeric. Trees with a height lower than this threshold (in meters)
#'   are flagged.
#' @param seeds A LAS object containing tree seeds, typically generated with
#'   \link{find_seeds}. If missing, an internal routine estimates seed positions
#'   from the lowest point of each tree.
#' @param buffer Numeric value (in meters) used to shrink the convex hull before
#'   filtering trees. Default is `-5`, which excludes trees within 5 meters of the
#'   boundary. Can also be an `sf` or `sfc` object for custom clipping. It must be a simple polygon,
#'   i.e. a single `POLYGON` without holes (multipolygons are not supported), with the same CRS as `las`.
#'
#' @return A modified LAS object in which the points of flagged trees have their `UserData`
#'   attribute set to a non-zero value. `treeID` is not modified.
#'
#' @seealso \link{qsf}, \link{find_seeds}
#'
#' @name flagging
#' @md
#' @export
flag_small_trees = function(las, max_height = 2)
{
  if (!"treeID" %in% names(las))   stop("Input point cloud must have an attribute 'treeID'")
  if (!"hag" %in% names(las))      stop("Input point cloud must have an attribute 'hag'")

  # C++ writes UserData in place: allocate a fresh copy to not modify the input
  if (!"UserData" %in% names(las)) las@data$UserData <- ARBORTREE
  las@data$UserData <- data.table::copy(las@data$UserData)
  flag_small_trees_cpp(las@data, max_height)
  las
}

keep_small_trees = function(las, max_height = 2)
{
  UserData <- treeID <- hag <- hag_max <- hag_min <- NULL

  attributes <- names(las)
  stopifnot("treeID" %in% attributes)
  stopifnot("hag" %in% attributes)

  ans   <- las@data[!is.na(treeID), list(hag_max = max(hag), hag_min = min(hag)), by = treeID]
  ans   <- ans[hag_max < max_height]

  trees <- lidR::filter_poi(las, treeID %in% ans$treeID & !is.na(treeID))

  trees
}


#' @rdname flagging
#' @export
#' @importFrom data.table :=
flag_buffer = function(las, seeds, buffer = -5)
{
  is_polygon <- inherits(buffer, "sf") || inherits(buffer, "sfc")

  if (is.numeric(buffer))
  {
    if (buffer > 0)
      stop("'buffer' must be negative.")
  }
  else if (is_polygon)
  {
    geom <- sf::st_geometry(buffer)
    if (length(geom) != 1L || sf::st_geometry_type(geom) != "POLYGON" || length(geom[[1]]) != 1L)
      stop("'buffer' must be a simple polygon: a single POLYGON without holes.")
  }
  else
  {
    stop("'buffer' must be a number or an sf polygon.")
  }

  if (!"treeID" %in% names(las))
  {
    stop("Input point cloud must have an attribute 'treeID'")
  }

  if (!missing(seeds))
  {
    if (!"treeID" %in% names(seeds))
    {
      stop("Input seeds must have an attribute 'treeID'")
    }
  }

  if (is_polygon && sf::st_crs(geom) != sf::st_crs(las))
    stop("'buffer' and 'las' must have the same CRS.")

  # C++ writes UserData in place: allocate a fresh copy to not modify the input
  if (!"UserData" %in% names(las)) las@data$UserData <- ARBORTREE
  las@data$UserData <- data.table::copy(las@data$UserData)

  see <- if (missing(seeds)) NULL else seeds@data

  if (is_polygon)
  {
    xy <- sf::st_coordinates(geom)
    flag_buffer_polygon_cpp(las@data, see, xy[, "X"], xy[, "Y"])
  }
  else
  {
    flag_buffer_cpp(las@data, see, buffer)
  }

  las
}

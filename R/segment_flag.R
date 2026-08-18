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
#' filter data on the fly. See the [Arbor book](<placeholder>)
#'
#' `flag_small_trees()` flags trees below a given height threshold, which are often
#' poorly segmented or correspond to understory vegetation.\cr\cr
#' `flag_buffer()` flags trees located near the edges of a point cloud by clipping
#' trees whose seed points fall outside a buffered polygon.\cr\cr
#' See the [Arbor book](<placeholder>) for more details.
#'
#' @param las A LAS object from lidR containing segmented trees.
#' @param max_height Numeric. Trees with a height lower than this threshold (in meters)
#'   are flagged.
#' @param seeds A LAS object containing tree seeds, typically generated with
#'   \link{find_seeds}. If missing, an internal routine estimates seed positions
#'   from the lowest point of each tree.
#' @param buffer Numeric value (in meters) used to shrink the convex hull before
#'   filtering trees. Default is `-5`, which excludes trees within 5 meters of the
#'   boundary. Can also be an `sf` polygon object for custom clipping.
#'
#' @return A modified LAS object with flagged trees assigned negative `treeID` values.
#'
#' @seealso \link{qsf}, \link{find_seeds}
#'
#' @name flagging
#' @md
#' @export
flag_small_trees = function(las, max_height = 2)
{
  UserData <- treeID <- hag <- hag_max <- hag_min <- NULL

  attributes <- names(las)
  if (!"treeID" %in% names(las))   stop("Input point cloud must have an attribute 'treeID'")
  if (!"hag" %in% names(las))      stop("Input point cloud must have an attribute 'hag'")

  ans <- las@data[treeID > 0, list(hag_max = max(hag), hag_min = min(hag)), by = treeID]
  ans <- ans[hag_max > max_height & hag_min < max_height]

  # Keep trees whose seeds are inside
  if (!"UserData" %in% names(las)) las@data$UserData = ARBORTREE
  las@data$UserData = data.table::copy(las@data$UserData)
  las@data[UserData == ARBORUNDERSTORY, UserData := ARBORTREE] # revert previous
  las@data[!treeID %in% ans$treeID & UserData != ARBORLOW, UserData := ARBORUNDERSTORY]
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
  if (buffer > 0) stop("'buffer' must be negative.")

  UserData <- NULL

  if (!"treeID" %in% names(las))   stop("Input point cloud must have an attribute 'treeID'")
  if (!missing(seeds))
  {
    if (!"treeID" %in% names(seeds)) stop("Input seeds must have an attribute 'treeID'")
  }

  # Avoid NOTES
  X <- Y <- Z <- treeID <- NULL

  # Compute roots per tree (drop unsegmented)
  if (missing(seeds))
  {
    # Make root robust to small groups
    root <- function(x, y, z)
    {
      i <- order(z)
      n <- min(100L, length(i))            # Use available points (<=100)
      if (n < 10L) {                       # Guard against tiny groups
        return(list(X = NA_real_, Y = NA_real_, Z = NA_real_))
      }
      x <- x[i][seq_len(n)]
      y <- y[i][seq_len(n)]
      z <- z[i][seq_len(n)]
      x <- mean(x)
      y <- mean(y)
      z <- mean(z)
      return(list(X = x, Y = y, Z = z))
    }

    roots <- las@data[treeID > 0, root(X, Y, Z), by = treeID]
    roots <- roots[is.finite(X) & is.finite(Y)]     # Drop invalid coords
    if (!nrow(roots)) return(las[0])                # Early exit if nothing valid
  }
  else
  {
    # Make root robust to small groups
    root <- function(x, y, z)
    {
      x <- mean(x)
      y <- mean(y)
      z <- mean(z)
      return(list(X = x, Y = y, Z = z))
    }

    roots <- seeds@data[treeID > 0, root(X, Y, Z), by = treeID]
    roots <- roots[is.finite(X) & is.finite(Y)]     # Drop invalid coords
    if (!nrow(roots)) return(las[0])                # Early exit if nothing valid
  }

  # Build buffered convex hull
  if (is.numeric(buffer))
  {
    bb <- suppressWarnings(sf::st_convex_hull(las))
    bb <- suppressWarnings(sf::st_buffer(bb, dist = buffer))
    if (any(sf::st_is_empty(bb))) return(las[0])    # Early exit if buffer kills hull
  }

  if (inherits(buffer, "sf") || inherits(buffer, "sfc"))
  {
    if (!any(sf::st_geometry_type(buffer) %in% c("POLYGON", "MULTIPOLYGON")))
      stop("buffer must be a polygon or multipolygon")

    bb <- buffer
  }

  # Convert roots to sf (safe: no NA)
  seeds <- sf::st_as_sf(roots, coords = c("X", "Y"), crs = sf::st_crs(las))
  idx <- sf::st_intersects(seeds, bb)
  valid_seeds <- seeds[lengths(idx) > 0, ]
  if (!nrow(valid_seeds)) return(las[0])

  # Keep trees whose seeds are inside
  if (!"UserData" %in% names(las)) las@data$UserData = ARBORTREE
  las@data$UserData = data.table::copy(las@data$UserData)
  las@data[UserData == ARBORBUFFER, UserData := ARBORTREE] # revert previous
  las@data[!treeID %in% valid_seeds$treeID & UserData != ARBORLOW, UserData := ARBORBUFFER]
  return(las)
}


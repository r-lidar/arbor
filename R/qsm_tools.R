# @file qsm_tools.R
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

#' Simple tools for QSM
#'
#' Simple tools for QSM.
#'
#' @param qs A QSM or QSF
#' @param qsm A QSM
#' @export
#' @rdname qsm_tools
#' @md
qsm_volume <- function(qs)
{
  UseMethod("qsm_volume")
}

#' @export
qsm_volume.qsm = function(qs)
{
  l = sqrt((qs$endX-qs$startX)^2+(qs$endY-qs$startY)^2+(qs$endZ-qs$startZ)^2)
  qs$volume = pi*qs$radius^2*l
  V = sum(qs$volume)
  return(V)
}

#' @export
qsm_volume.qsf <- function(qs)
{
  vapply(qs, qsm_volume, numeric(1))
}

#' @rdname qsm_tools
#' @export
qsm_height <- function(qs)
{
  UseMethod("qsm_height")
}

#' @export
qsm_height.qsm <- function(qs)
{
  if (nrow(qs) == 0) return(NA_real_)
  zmin <- min(c(qs$startZ, qs$endZ), na.rm = TRUE)
  zmax <- max(c(qs$startZ, qs$endZ), na.rm = TRUE)
  zmax - zmin
}

#' @export
qsm_height.qsf <- function(qs)
{
  vapply(qs, qsm_height, numeric(1))
}

#' @rdname qsm_tools
#' @param short bool. Shorter message: only the human readable label is returned
#' (the optional leading `[<code>]` tag, see \link{qsf_log}, is stripped as well).
#' @export
qsm_message = function(qsm, short = FALSE)
{
  msg = attr(qsm, "message")
  if (short) msg = vapply(msg, function(m) .qsf_parse_message(m)$label, character(1), USE.NAMES = FALSE)
  msg
}

filter_tree = function(tree)
{
  foliage <- NULL
  attributes = names(tree)

  if ("foliage" %in% attributes)
  {
    tree = lidR::filter_poi(tree, foliage == FALSE)
  }

  if ("treeID" %in% attributes)
  {
    if (length(unique(tree$treeID)) != 1)
      stop("The point cloud must contain a single tree", call. = FALSE)
  }

  return(tree)
}

#' Add Synthetic Ground Points to a Single-Tree Point Cloud
#'
#' Adds a set of synthetic ground points below a single-tree point cloud. Useful to perform
#' arbor operations on already externally isolated trees
#'
#' @param las A `LAS` object containing a single-tree point cloud.
#' @param n integer. Number of ground points
#'
#' @return A `LAS` object containing the original points plus n synthetic
#' ground points.
#'
#' @export
#' @md
add_single_tree_ground = function(las, n = 1000)
{
  z  <- min(las$Z)
  bb <- lidR::st_bbox(las)

  xg <- stats::runif(n, bb[1] - 1, bb[3] + 1)
  yg <- stats::runif(n, bb[2] - 1, bb[4] + 1)

  lidR::quantize(xg, 0.001, las@header[["X offset"]])
  lidR::quantize(yg, 0.001, las@header[["Y offset"]])

  # Create a data.frame with same columns as input LAS
  template <- las@data[rep(1, n), , drop = FALSE]

  # Reset all values
  for (col in names(template))
  {
    if (is.integer(template[[col]]))
      template[[col]] <- NA_integer_
    else if (is.numeric(template[[col]]))
      template[[col]] <- NA_real_
    else if (is.logical(template[[col]]))
      template[[col]] <- NA
    else
      template[[col]] <- NA
  }

  # Fill mandatory coordinates
  template$X <- xg
  template$Y <- yg
  template$Z <- z

  # Set common ground defaults if present
  if ("Classification" %in% names(template))
    template$Classification <- 2L  # ASPRS ground

  if ("ReturnNumber" %in% names(template))
    template$ReturnNumber <- 1L

  if ("NumberOfReturns" %in% names(template))
    template$NumberOfReturns <- 1L

  gnd <- suppressWarnings(lidR::LAS(template, header = las@header))
  suppressWarnings(rbind(las, gnd))
}






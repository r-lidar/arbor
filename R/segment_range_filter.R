# @file segment_range_filter.R
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

#' Trajectory-based range filtering utilities
#'
#' These functions implement a simple workflow to:
#' \enumerate{
#'   \item read trajectory file,
#'   \item compute point-to-sensor range using the trajectory,
#'   \item filter points based on range and height above ground.
#' }
#' They are designed to be used together, but can also be called independently.
#'
#' @name trajectory_range
#' @rdname trajectory_range
NULL


#' Read a trajectory file
#'
#' `read_trajectory`: Reads a trajectory file and converts it to an \code{sf} point object.
#' The first column is assumed to be the GPS time and is renamed to \code{gpstime}.
#' The file must contain columns named \code{x}, \code{y}, and \code{z}.
#'
#' @param file Character. Path to the ASCII trajectory file
#' @param dt Numeric. Keep a record every \code{dt} seconds. Default 0.5.
#' @return An \code{sf} object with point geometry and a \code{gpstime} attribute.
#'
#' @md
#' @export
#' @rdname trajectory_range
read_trajectory = function(file, dt = 0.5)
{
  traj <- data.table::fread(file)
  if (ncol(traj) < 4) stop("Trajectory file must contain at least gpstime, x, y and z columns")
  names(traj)[1] <- "gpstime"
  traj = sf::st_as_sf(traj, coords = c("x", "y", "z"))

  if (dt > 0)
    traj = traj[ !duplicated(floor(traj$gpstime / dt)), ]

  return(traj)
}


#' Add point-to-sensor range attribute
#'
#' `add_range`: Computes the distance between each point and the sensor position
#' derived from a trajectory, and stores it as a new LAS attribute named
#' \code{range}.
#'
#' @param las A \code{LAS} object.
#' @param traj An \code{sf} trajectory object as returned by \link{read_trajectory}.
#'
#' @return A \code{LAS} object with an additional attribute \code{range}.
#'
#' @md
#' @export
#' @rdname trajectory_range
add_range = function(las, traj)
{
  range_data <- lidR::get_range(las, traj)
  lidR::add_lasattribute(las, range_data, "range", "distance to sensor")
}


#' Filter points by range and height above ground
#'
#' `filter_range`: Filters a point cloud by removing points that are too far to the sensor
#' and below a given height above ground. For example, if \code{distance = 15}
#' and \code{hag_max = 4}, points that are beyond 15 m of the sensor \emph{and}
#' lower than 4 m above ground are removed. Points farther than 15 m from the
#' sensor but higher than 4 m above ground, are retained. This prevents the
#' removal of canopy point.
#'
#' @param distance Numeric. Maximum allowed distance to the sensor.
#' @param hag_max Numeric. Maximum allowed height above ground.
#'
#' @return A filtered \code{LAS} object.
#'
#' @export
#' @md
#' @rdname trajectory_range
filter_range <- function(las, distance = 15, hag_max = 4)
{
  range <- hag <- NULL
  if (!"hag" %in% names(las)) stop("The point cloud must have an 'hag' attribute")
  n_before <- lidR::npoints(las)
  las <- lidR::filter_poi(las, (hag >= hag_max) | (range <= distance))
  n_after <- lidR::npoints(las)
  cat("Removed", round(100 * (n_before - n_after) / n_before, 2), "% of points\n")
  las
}

# @file segment_ground.R
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

#' Ground Segmentation
#'
#' Classifies ground points and computes height above ground. See
#' the [Arbor book](<placeholder>) for mode details.
#'
#' @param las A `LAS` object (from the `lidR` package) containing the point cloud data
#' @param params list See \link{parameters}.
#'
#' @md
#' @export
#' @seealso \link{find_seeds}, \link{segment_semantic}, \link{segment_instance}
segment_ground = function(las, params = arbor_parameters_default)
{
  if (!"Classification" %in% names(las)) las@data$Classification = 0L
  las <- lidR::add_lasattribute(las, 0, "hag", "Height Above Ground")
  segment_ground_cpp(las@data, params)
  return(las)
}

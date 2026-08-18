# @file segment_decimation.R
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

#' Homogenization of the Point Cloud
#'
#' Homogenization of the point cloud using a hybrid approach that includes
#' Barycentric Voxel Decimation and reinjection of random points. See
#' the [Arbor book](https://r-lidar.github.io/arbor_book/) for mode details.
#'
#' @param las LAS object from lidR
#' @param res Voxel resolution.
#' @md
#' @export
hybrid_homogeneization = function(las, res = 0.02)
{
  stop_if_not_tls(las)
  decimated <- TRUE
  keep <- C_homogeneization(las@data, res)
  las <- las[keep]
  free(keep)
  return(las)
}

barycentric_decimation = function(las, res)
{
  decimated <- TRUE
  keep <- C_homogeneization(las@data, res, hybrid = FALSE)
  las <- las[keep]
  free(keep)
  return(las)
}

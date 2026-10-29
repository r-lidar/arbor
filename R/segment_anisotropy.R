# @file segment_anisotropy.R
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

#' Compute Local Wood Likelihood
#'
#' This function calculates the local anisotropy of each point in a LAS object using
#' a k-nearest neighbors (KNN) approach. If the point cloud has an attribute 'hag', the points at or
#' below the cut height are not computed and get 0: the semantic segmentation classifies them foliage
#' anyway. See the [Arbor book](https://r-lidar.github.io/arbor_book/) for mode details.
#'
#' @param las A LAS object from lidR.
#' @param params list See \link{parameters}.
#' @export
wood_likelihood = function(las, params = arbor_parameters_default)
{
  stop_if_not_tls(las)

  # The wood likelihood is written in place by C++: R must allocate the attribute 'pwood' first
  las@data$pwood <- 0
  las <- lidR::add_lasattribute_manual(las, name = "pwood", desc = "wood likelyhood", type = "float")

  wood_likelihood_cpp(las@data, params)
  return(las)
}

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
#' a k-nearest neighbors (KNN) approach. See the [Arbor book](https://r-lidar.github.io/arbor_book/) for mode details.
#'
#' @param las A LAS object from lidR.
#' @param params list See \link{parameters}.
#' @export
wood_likelihood = function(las, params = arbor_parameters_default)
{
  stop_if_not_tls(las)
  k <- params$woodlikelihood$k
  anisotropy <- C_anisotropy(las@data, k)
  las <- lidR::add_lasattribute_manual(las, anisotropy, "pwood", "wood likelyhood", "float")
  free(anisotropy)
  return(las)
}

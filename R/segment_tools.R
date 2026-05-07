# @file segment_tools.R
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

evaluate_penalty = function(params)
{
  penalty = params$path_finder$angle_penalty(0:180)
  if (length(penalty) != 181) stop("Invalid penalty function")
  params$path_finder$penalty = penalty
  params
}

connected_components = function(las, res, min_pts, name = "clusterID", connectivity = 6)
{
  .N <- N <- clusterID <- NULL

  allowed <- c(6L, 18L, 26L)
  if (!connectivity %in% allowed) stop(sprintf("Invalid connectivity: %d. Allowed values are %s", connectivity, paste(allowed, collapse = ", ")))

  u = C_connected_component(las@data, res, connectivity)
  las = lidR::add_lasattribute(las, u, name, "connected component ID")
  grp = las@data[, .N, by = clusterID]
  grp = grp[N < min_pts]
  invalid = las@data[[name]] %in% grp$clusterID
  las@data[[name]][invalid] = 0L
  return(las)
}

sor = function(las, k, m)
{
  noise = C_sor(las@data, k, m)

  if ("Classification" %in% names(las))
  {
    new_classes <- las@data[["Classification"]]
    new_classes[new_classes == lidR::LASNOISE] <- lidR::LASUNCLASSIFIED
  }
  else
  {
    new_classes <- rep(lidR::LASUNCLASSIFIED, lidR::npoints(las))
  }

  new_classes[noise] <- lidR::LASNOISE
  las@data[["Classification"]] <- new_classes
  return(las)
}


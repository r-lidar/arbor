# @file segment_transfert_attributes.R
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

#' Transfer attributes between point clouds
#'
#' Transfers one or more attributes from a decimated (sparser) point cloud
#' to a denser point cloud using a 1-nearest neighbor (1-NN) approach.
#' Each point in the target point cloud receives the attribute value of
#' its nearest neighbor in the source point cloud.
#'
#' @param from A `LAS` object representing the source (typically decimated) point cloud
#'   containing the attributes to transfer.
#' @param to A `LAS` object representing the target (typically denser) point cloud
#'   that will receive the transferred attributes.
#' @param attributes A character vector specifying the names of the attributes
#'   to transfer from `from` to `to`.
#'
#' @details
#' The function performs a 1-NN search between the `from` and `to` point clouds
#' using `lidR::knnx()`. For each attribute listed in `attributes`, the values are
#' assigned from the nearest point in `from`. If available, the function also
#' propagates the corresponding attribute description from the Extra Bytes
#' metadata (`VLR$Extra_Bytes`).
#'
#' @return A `LAS` object identical to `to`, but with the requested attributes added.
#'
#' @seealso [lidR::knnx()], [lidR::add_lasattribute()]
#'
#' @examples
#' \dontrun{
#' las <- transfer_attributes(las_decimated, las, attributes = "treeID", "foliage")
#' }
#'
#' @export
transfer_attributes = function(from, to, attributes)
{
  # Sanity check
  stopifnot(all(attributes %in% names(from@data)))

  # Compute 1-NN mapping (index + distance)
  nn   <- lidR::knnx(from, to, k = 1)
  idx  <- nn[[1]]
  dist <- nn[[2]]

  # For each requested attribute, transfer it
  for (attr in attributes)
  {
    eba <- from@header@VLR$Extra_Bytes[["Extra Bytes Description"]][[attr]]
    values <- from@data[[attr]][idx]
    to@data[[attr]] <- values

    if (!is.null(eba))
      to <- lidR::add_lasattribute(to,  name = attr,  desc = eba$description)
  }

  return(to)
}

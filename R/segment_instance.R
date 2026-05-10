# @file segment_instance.R
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

#' Individual Tree Instance Segmentation
#'
#' This function performs individual tree instance segmentation by constructing a point-based network
#' and computing the least-cost path for each point to predefined seed points. The goal is to assign
#' each point in the point cloud a `treeID` corresponding to its nearest seed, based on the least-cost
#' path computed through the network. The process incorporates both spatial and heuristic constraints
#' to segment trees effectively. See the [Arbor book](<placeholder>) for mode details.
#'
#' @param las A `LAS` object (from the `lidR` package) containing the point cloud data. It must have
#'   an attribute `foliage`, which can be generated using prior semantic segmentation.
#' @param seeds An `sf` object containing points that represent tree seeds. Each point must have
#'   xyz coordinates and a `treeID`. Seeds can be generated using \link{find_seeds}.
#' @param params list See \link{parameters}.
#'
#' @return A `LAS` object with an additional attribute `treeID`, indicating the ID of the tree to
#'   which each point belongs.
#'
#' @md
#' @export
#' @rdname segment_instance
#' @name segment_instance
#' @importFrom Rcpp sourceCpp
#' @seealso \link{find_seeds}, \link{segment_semantic}
segment_instance = function(las, seeds, params)
{
  stop_if_not_tls(las)
  las@data$treeID = NA_integer_
  segment_instance_cpp(las@data, seeds@data, params)
  las <- lidR::add_lasattribute(las, name = "treeID", desc = "Unique ID per tree")
  return(las)
}


# segment_instance_r = function(las, seeds, params)
# {
#   logger("Instance segmentation start")
#
#   # The point cloud must have hag, and foliage computed
#   attributes = names(las)
#   stopifnot("hag" %in% attributes)
#   stopifnot("foliage" %in% attributes)
#   if (!"pointID" %in% names(las)) las@data$pointID = 1:lidR::npoints(las)
#
#   logger("Decimating the point cloud... (1/4)")
#
#   res        <- params$path_finder$decimation
#   core       <- hybrid_homogeneization(las, res)
#
#   num_points <- lidR::npoints(core)
#   num_trees  <- nrow(seeds)
#
#   cat(sprintf("raw %d, dec %d, seeds %d\n", lidR::npoints(las), num_points, num_trees));
#
#   # Plot for debugging
#   if (FALSE)
#   {
#     x <- plot(core)
#     plot(seeds, add = x, pal = "green", size = 6)
#   }
#
#   logger("Constructing the graph object (2/4)")
#
#   params <- evaluate_penalty(params)
#   graph  <- build_instance_graph(core@data, seeds@data, params);
#
#   logger("Pathfinder (3/4)")
#
#   seeds_ids  <- (num_points):(num_points+num_trees-1)
#   treeID     <- find_closest_node(graph, seeds_ids)
#
#   trueTreeID <- treeID[1:lidR::npoints(core)]
#   trueTreeID <- trueTreeID - min(seeds_ids) + 1
#   ID         <- seeds$treeID[trueTreeID]
#   core       <- lidR::add_lasattribute(core, ID, name = "treeID", desc = "tree ID")
#
#   if (FALSE)
#   {
#     x = plot(core, color = "treeID", size =2)
#     plot(seeds, add = x, pal = "red", size = 6)
#   }
#
#   logger("Assigning tree IDs to the dense point cloud (4/4)")
#
#   las <- transfer_attributes(core, las, "treeID")
#
#   logger("Instance segmentation completed")
#
#   return(las)
# }

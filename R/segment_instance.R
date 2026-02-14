#' Individual Tree Instance Segmentation
#'
#' This function performs individual tree instance segmentation by constructing a point-based network
#' and computing the least-cost path for each point to predefined seed points. The goal is to assign
#' each point in the point cloud a `treeID` corresponding to its nearest seed, based on the least-cost
#' path computed through the network. The process incorporates both spatial and heuristic constraints
#' to segment trees effectively.\cr\cr
#' The segmentation algorithm operates as follows:
#' \enumerate{
#' \item **Network Construction**: A k-nearest neighbor (k-NN) graph is built where points are connected
#' to their `k` closest neighbors. The cost to move in the graph between two points is explained in
#' next section and is calculated based on criteria such as euclidean distance, foliage type, and
#' movement direction
#' \item **Pathfinding**: Using the constructed graph, the function computes the least-cost path from
#' every point to each seed.
#' \item **Tree ID Assignment**: Each point is assigned the `treeID` of the seed with the least-cost
#' path.
#' }
#' Key Features of the Cost Function:
#' \itemize{
#' \item **Proximity Penalty**: Close points are favored with lower movement costs by using the cube of
#' the euclidean distance
#' \item **Foliage Type Penalty**: Foliage to foliage movement or foliage to wood movement incur
#' substantial additional costs.
#' \item **Directional Penalty**: Movements against gravity (upward) are penalized more to prioritize
#' paths leading downward to seeds.
#' }
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
#' @details
#' Internal Processing Steps:
#' 1. **Preprocessing**:
#'    - Decimate the point cloud based on the specified resolution (`res`).
#'    - Transform and scale coordinates for network construction and pathfinding.
#' 2. **Network Construction**:
#'    - Points are connected to their `k` nearest neighbors with a distance threshold (`max_gap`).
#'    - The cost on each edge is the cube of the euclidean distance
#'    - Additional cost penalties are applied based on movement direction and wood/foliage transitions.
#' 3. **Pathfinding**:
#'    - A graph is built combining point and seed networks.
#'    - The shortest path (least-cost) from each point to all seeds is calculated.
#' 4. **Tree Assignment**:
#'    - Points are assigned the `treeID` of the nearest seed based on path costs.
#'    - IDs are propagated back to the original dense point cloud.
#'
#' @md
#' @export
#' @rdname segment_instance
#' @name segment_instance
#' @importFrom Rcpp sourceCpp
#' @seealso \link{find_seeds}, \link{segment_semantic}
segment_instance = function(las, seeds, params)
{
  params <- evaluate_penalty(params)
  treeID <- segment_instance_cpp(las@data, seeds@data, params)
  las    <- lidR::add_lasattribute(las, treeID, name = "treeID", desc = "Unique ID per tree")
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

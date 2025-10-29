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
#' @importFrom Rcpp sourceCpp
#' @seealso \link{find_seeds}, \link{segment_foliage}
segment_vegetation = function(las, seeds, params)
{
  #res = 0.05; k = 5; max_gap = 0.5; z_factor = 0.8
  z_factor <- params$instance$z_scale

  . <- X <- Y <- Z <- foliage <- pointID <- NULL

  cost_factors = list(
    wood2wood = 0.1,
    leaf2leaf = 20,
    wood2leaf = 1000)

  # The point cloud must have hag, and foliage computed
  attributes = names(las)
  stopifnot("hag" %in% attributes)
  stopifnot("foliage" %in% attributes)

  ti = tic()

  if (!"pointID" %in% names(las)) las@data$pointID = 1:lidR::npoints(las)

  cat("Decimating the point cloud... (1/6)\n") ; t0 = tic()

  # We will apply path finding to a decimated point cloud
  ans        <- decimate_translate(las, params)
  dec        <- ans$dec
  dec$Z      <- dec$Z * z_factor
  dec@data   <- dec@data[, .(X,Y,Z, foliage, pointID)]
  num_points <- lidR::npoints(dec)

  # Global translation to origin for computation stability
  # (translation already applied to dec)
  x_translation <- ans$x_translation
  y_translation <- ans$y_translation

  # Get the coordinates of the seeds
  seeds@data$X <- seeds@data$X - x_translation
  seeds@data$Y <- seeds@data$Y - y_translation
  seeds@data$Z <- seeds@data$Z * z_factor
  num_trees    <- nrow(seeds)

  # Compute the single master seed that rule them all. the real spatial position
  # does not matter since the cost is 0 but having a position is better for rendering
  master_seed <- make_master_seed(dec)

  # Plot for debugging
  if (FALSE)
  {
    x <- plot(dec)
    plot(seeds, add = x, pal = "green", size = 6)
    plot(master_seed, add = x, pal = "white", size = 8)
  }


  toc(t0)
  cat("Constructing the graph object... (Step 4/6)\n") ; t0 = tic()

  # Each point is connected to its knn. The connection is bidirectional. The cost of the connection
  # is based on the euclidean distance with some variation in order to specifically follow the wood
  # and not foliage. The cost is not the same in both directions
  k <- params$path_finder$k_neighborhood_connectivity
  max_gap <- params$path_finder$max_gap

  graph <- build_instance_graph(dec@data, seeds@data, master_seed@data, k = k, max_gap = max_gap);

  toc(t0)
  cat("Pathfinder... (Step 5/6)\n") ; t0 = tic()

  seeds_ids <- (num_points):(num_points+num_trees-1)
  ans = find_closest_ground(graph, seeds_ids)

  treeID     <- ans$closest_ground
  trueTreeID <- treeID[1:lidR::npoints(dec)]
  trueTreeID <- trueTreeID - min(seeds_ids) + 1 # because there is an index error somewhere
  ID         <- seeds$treeID[trueTreeID]
  dec <- lidR::add_lasattribute(dec, ID, name = "treeID", desc = "tree ID")

  toc(t0)
  cat("Assigning tree IDs to the dense point cloud... (Step 6/6)\n") ; t0 = tic()

  dec$X <- dec$X + x_translation
  dec$Y <- dec$Y + y_translation
  dec$Z <- dec$Z / z_factor

  las <- transfer_attributes(dec, las, "treeID")

  toc(ti, space = "")
  gc()

  return(las)
}

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
#' @param res Numeric. The resolution for decimating the point cloud. Smaller resolutions retain more
#'   points but increase computational cost. Ideally keeping one point every 5 cm would be optimal.
#'   In practice it increases the computation times. Every 10 cm works well actually and compute fast.
#'   The default is 0.08 meters (8 cm), a trade-off between 5 and 10 cm.
#' @param max_gap Numeric. Maximum allowed distance (in meters) between points for them to be connected
#'   in the network. Points separated by larger gaps are disconnected, effectively creating barriers.
#'   Default is 0.2 meters.
#' @param k Integer. Number of nearest neighbors used to build the k-NN graph for the network. Default
#'  is 10.
#' @param z_factor Numeric. Scaling factor for the Z-axis. The point cloud is compressed vertically to
#'   reduce the effect of large vertical gaps in the canopy while preserving horizontal distances.
#'   Default is 0.8.
#' @param ... unused. Serves only to separate easy parameters to complex ones that should not need to
#' be modified.
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
segment_vegetation = function(las, seeds, ..., res = 0.05, k = 5, max_gap = 0.5, z_factor = 0.8)
{
  #res = 0.05; k = 5; max_gap = 0.5; z_factor = 0.8

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
  ans <- decimate_translate(las, res, z_factor)
  dec <- ans$dec
  dec@data <- dec@data[, .(X,Y,Z, foliage, pointID)]
  num_points <- lidR::npoints(dec)

  # Global translation to origin for computation stability
  # (translation already applied to dec)
  x_translation <- ans$x_translation
  y_translation <- ans$y_translation

  # Get the coordinates of the seeds
  seeds@data$X <- seeds@data$X - x_translation
  seeds@data$Y <- seeds@data$Y - y_translation
  seeds@data$Z <- seeds@data$Z * z_factor
  num_trees <- nrow(seeds)

  # Compute the single master seed that rule them all. the real spatial position
  # does not matter since the cost is 0 but having a position is better for rendering
  master_seed <- make_master_seed(dec)

  # Plot for debugging
  if (FALSE)
  {
    x = plot(dec)
    plot(seeds, add = x, pal = "green", size = 6)
    plot(master_seed, add = x, pal = "white", size = 8)
  }

  toc(t0)
  cat("Identifying neighbors... (2/6)\n") ; t0 = tic()

  # Each point is connected to its knn. The connection is bidirectional. The cost of the connection
  # is based on the euclidean distance with some variation in order to specifically follow the wood
  # and not foliage. The cost is not the same in both directions

  point_network <- compute_point_network(dec, k = k, max_gap = max_gap, wood_mask = dec$foliage, cost_factors = cost_factors)
  points_ids = 1:num_points

  toc(t0)
  cat("Building target connectivity... (3/6)\n") ; t0 = tic()

  # Each seed is connected to the knn in the point cloud. The connection is undirectional.
  # it is possible to move from the seed to the scene not the opposite. The cost is the distance
  seed_network  <- compute_network(dec, seeds, k = k)
  seed_network$from <- seed_network$from + num_points
  seeds_ids <- (num_points+1):(num_points+num_trees)

  master_seed_network  <- compute_network(seeds, master_seed, k = num_trees)
  master_seed_network$from <- master_seed_network$from + max(seeds_ids)
  master_seed_network$to <- master_seed_network$to + min(seeds_ids)
  master_seed_network$cost <- 0.001
  master_seed_id <- master_seed_network[1,1]

  toc(t0)
  cat("Constructing the graph object... (Step 4/6)\n") ; t0 = tic()


  combined_network <- rbind(point_network, seed_network, master_seed_network)
  free(point_network, seed_network, master_seed_network)

  graph <- build_graph(combined_network)
  cache <- compute_distances(graph, master_seed_id)

  toc(t0)
  cat("Pathfinder... (Step 5/6)\n") ; t0 = tic()

  from <- rep(master_seed_id, length(points_ids))
  to   <- points_ids

  # For loop by chunk to reduce memory usage and have an estimated progression
  chunk_size <- 50000
  chunks <- split(to, ceiling(seq_along(to) / chunk_size))
  pb <- utils::txtProgressBar(min = 0, max = length(chunks), style = 3, width = 50)

  treeID <- rep(NA_integer_, lidR::npoints(dec))

  for (i in seq_along(chunks))
  {
    current_to <- chunks[[i]]

    path = findPaths(graph, cache, from[seq_along(current_to)], current_to)
    path = path$paths
    path <- lapply(path, function(x) x[2])
    tree_id_vector <- unlist(path)
    treeID[current_to] = tree_id_vector
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  free(combined_network)

  trueTreeID = treeID - min(seeds_ids)  +1 #because there is an index error somewhere
  ID = seeds$treeID[trueTreeID]
  dec <- lidR::add_lasattribute(dec, ID, name = "treeID", desc = "tree ID")

  toc(t0)
  cat("Assigning tree IDs to the dense point cloud... (Step 6/6)\n") ; t0 = tic()

  dec$X <- dec$X + x_translation
  dec$Y <- dec$Y + y_translation
  dec$Z <- dec$Z / z_factor

  las <- expand_treeid_to_neighbors(las, dec, z_factor = z_factor)

  toc(ti, space = "")
  gc()

  return(las)
}

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
segment_vegetation = function(las, seeds, ..., res = 0.08, k = 10, max_gap = 0.5, z_factor = 0.8)
{
  # Other parameters are hard coded
  upward_cost_factor = 100
  wood2wood_cost_factor = 0.1
  leaf2leaf_cost_factor = 20
  wood2leaf_cost_factor = 1000

  # The point cloud must have hag, anisotropy and foliage computed
  attributes = names(las)
  stopifnot("anisotropy" %in% attributes)
  stopifnot("foliage" %in% attributes)

  ti = tic()

  seed_coordinates <- sf::st_coordinates(seeds)

  seeds <- as.data.frame(seeds)
  seeds$X <- seed_coordinates[, "X"]
  seeds$Y <- seed_coordinates[, "Y"]
  seeds$Z <- seed_coordinates[, "Z"]

  x = mean(las$X)
  y = mean(las$Y)
  z = mean(las$Z)-100
  master_seed = data.frame(X = x, Y = y, Z = z)
  lidR::quantize(master_seed[["X"]], 0.01, las@header$`X offset`)
  lidR::quantize(master_seed[["Y"]], 0.01, las@header$`Y offset`)
  lidR::quantize(master_seed[["Z"]], 0.01, las@header$`Z offset`)
  header = rlas::header_create(master_seed)
  master_seed = suppressWarnings(lidR::LAS(master_seed, header))

  cat("Decimating the point cloud... (Step 1/7)\n")

  t0 = tic()

  if (!"decimated" %in% names(las))
    dec = lidR::decimate_points(las, lidR::barycenter_per_voxel(res))
  else
    dec = lidR::filter_poi(las, decimated == TRUE)

  # Calculate translation to origin
  x_translation <- mean(las$X)
  y_translation <- mean(las$Y)

  dec@data$X <- dec@data$X - x_translation
  dec@data$Y <- dec@data$Y - y_translation
  dec@data$Z <- dec@data$Z * z_factor

  header = rlas::header_create(seeds)
  seeds = lidR::LAS(seeds[, c("X", "Y", "Z", "treeID")], header)
  seeds@data$X <- seeds@data$X - x_translation
  seeds@data$Y <- seeds@data$Y - y_translation
  seeds@data$Z <- seeds@data$Z * z_factor

  num_trees <- nrow(seeds)
  num_points <- nrow(dec)

  toc(t0)

  cat("Identifying neighbors... (Step 2/7)\n") ; t0 = tic()

  # Each point is connected to its knn. The connection is bidirectional. The cost of the connection
  # is based on the euclidean distance with some variation in order to specifically follow the wood
  # and not foliage. The cost is not the same in both directions

  point_network <- compute_network(dec, k = k)
  gc()

  # Gaps have an infinite cost
  gaps = point_network$cost > max_gap
  point_network$cost[gaps] = Inf
  free(gaps)

  # It is more expensive to move in large steps
  point_network$cost <- point_network$cost^2

  # It is more expensive to move upward  (downward actually because we are starting from the seeds)
  X1 = dec$X[point_network$from]
  X2 = dec$X[point_network$to]
  Y1 = dec$Y[point_network$from]
  Y2 = dec$Y[point_network$to]
  Z1 = dec$Z[point_network$from]
  Z2 = dec$Z[point_network$to]
  dx = X1-X2
  dy = Y1-Y2
  dz = Z1-Z2
  magnitude = sqrt(dx^2+dy^2+dz^2)
  dot_product <- -dz
  cos_theta = dot_product/magnitude
  angle_degree <- acos(cos_theta)*180/pi

  f = function(x)
  {
    y = exp(log(100)/100*x)
    y[x > 100] = 100
    y
  }

  #upward = Z1-Z2 > 0
  #point_network$cost[upward] = point_network$cost[upward]*upward_cost_factor
  point_network$cost = point_network$cost*f(angle_degree)
  free(dx, dy, dz, X1, X2, Y1, Y2, Z1, Z2, magnitude, dot_product, cos_theta, angle_degree)

  is_wood1 = !dec$foliage[point_network$from]
  is_wood2 = !dec$foliage[point_network$to]
  wood2wood = is_wood1 & is_wood2
  leaf2leaf = !is_wood1 & !is_wood2
  wood2leaf = is_wood1 & !is_wood2
  point_network$cost[wood2wood] = point_network$cost[wood2wood] * wood2wood_cost_factor
  point_network$cost[leaf2leaf] = point_network$cost[leaf2leaf] * leaf2leaf_cost_factor
  point_network$cost[wood2leaf] = point_network$cost[wood2leaf] * wood2leaf_cost_factor
  free(is_wood1, is_wood2, wood2wood, wood2leaf)

  points_ids = 1:num_points

  # Each seed is connected to the knn in the point cloud. The connection is undirectional.
  # it is possible to move from the seed to the scene not the opposite. The cost is the distance
  seed_network  <- compute_network(dec, seeds, k = k)
  seed_network$from <- seed_network$from + num_points
  seed_network$cost <- 0.001
  seeds_ids <- (num_points+1):(num_points+num_trees)

  master_seed_network  <- compute_network(seeds, master_seed, k = num_trees)
  master_seed_network$from <- master_seed_network$from + max(seeds_ids)
  master_seed_network$to <- master_seed_network$to + min(seeds_ids)
  master_seed_network$cost <- 0.001
  master_seed_id = master_seed_network[1,1]

  toc(t0)

  cat("Constructing the graph object... (Step 3/7)\n") ; t0 = tic()

  # Extract unique tree location IDs

  combined_network <- rbind(point_network, seed_network, master_seed_network)
  free(point_network, seed_network)

  #graph_object <- cppRouting::makegraph(combined_network, directed = TRUE)
  #free(combined_network)

  toc(t0)

  from = master_seed_id
  to = points_ids

  cat("Path finder... (Step 4/7)\n") ; t0 = tic()

  from = rep(master_seed_id, length(points_ids))
  to = points_ids

  t = tic()
  # For loop by chunk to reduce memory usage and have an estimated progression
  chunk_size <- 1000000
  chunks <- split(to, ceiling(seq_along(to) / chunk_size))
  pb <- utils::txtProgressBar(min = 0, max = length(chunks), style = 3, width = 50)

  treeID <- rep(NA_integer_, lidR::npoints(dec))

  for (i in seq_along(chunks))
  {
    current_to <- chunks[[i]]

    #path <- cppRouting::get_path_pair(graph_object, from = from[seq_along(current_to)], to = current_to)
    #path <- lapply(path, function(x) as.integer(x)[-1])
    #path <- unname(do.call(c, path))

    path = findPaths(combined_network, from[seq_along(current_to)], current_to)
    path = path$paths
    path <- lapply(path, function(x) x[2])
    tree_id_vector <- unlist(path)
    treeID[current_to] = tree_id_vector
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  toc(t)

  free(combined_network)

  # chunk_size <- 100000
  # chunks <- split(to, ceiling(seq_along(to) / chunk_size))
  # pb <- txtProgressBar(min = 0, max = length(chunks), style = 3, width = 50)
  # treeID <- rep(NA_integer_, npoints(dec))
  #
  # for (i in seq_along(chunks))
  # {
  #   current_to <- chunks[[i]]
  #   t = tic()
  #   distance_matrix <- cppRouting::get_distance_matrix(graph_object, from = seeds_ids, to = current_to)
  #   toc(t)
  #
  #   tree_id_vector = apply(distance_matrix, 2, which.min)
  #   ids = sapply(tree_id_vector, function(x) {  if (length(x) == 0) return(NA_integer_) else return(x) })
  #   ids = unname(unlist(ids))
  #   treeID[current_to] = seeds$treeID[ids]
  #
  #   setTxtProgressBar(pb, i)
  # }
  # close(pb)

  trueTreeID = treeID - min(seeds_ids) + 1 # +1 because there is an index error somewhere
  treeID = seeds$treeID[trueTreeID]
  dec <- lidR::add_lasattribute(dec, treeID, name = "treeID", desc = "tree ID")

  toc(t0)


  #plot(dec, color = "treeID", legend = TRUE)

  cat("Assigning tree IDs to the dense point cloud... (Step 6/7)\n") ; t0 = tic()

  dec$X <- dec$X + x_translation
  dec$Y <- dec$Y + y_translation
  dec$Z <- dec$Z / z_factor

  las = expand_treeid_to_neighbors(las, dec, z_factor = z_factor)

  toc(ti, space = "")

  gc()

  return(las)
}

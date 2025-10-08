#' Wood-Foliage Semantic Segmentation
#'
#' This function performs wood-foliage semantic segmentation by constructing a network of points
#' using k-nearest neighbors (KNN). A pathfinding algorithm navigates this network to identify the
#' tree skeleton by finding the least-cost path from every point in space to the ground. Each point
#' is classified as wood or foliage based on two main criteria: (1) Proximity to the detected skeleton and
#' (2) Anisotropy exceeding a given threshold. The function \link{compute_anisotropy} must be applied
#' first. The point cloud must have an attribute 'hag' (see \link[lidR:height_above_ground]{height_above_ground})
#' Finally, a connected component step removes small clusters incorrectly classified as wood,
#' since some foliage points may exhibit high anisotropy.\cr\cr
#' This function has multiple parameters. Those beyond `...` generally have a minor impact
#' and should only be adjusted in specific cases. If the scene contains regular trees, the
#' default values should be sufficient.
#'
#' @param las A LAS object from lidR.
#' @param dtm A SpatRaster object. Digital Terrain Model
#' @param res Resolution (in meters). The point cloud is initially downsampled to reduce
#'   computational cost. Ideally, retaining one point every 5 cm provides optimal results,
#'   but in practice, 10 cm works well and computes faster. The default is 8 cm as a
#'   trade-off between accuracy and speed.
#' @param max_gap Maximum allowable distance (in meters) between connected points in the
#'   network. Points farther apart than this threshold will not be linked. Default: 20 cm.
#' @param ... Unused. Serves only to distinguish common parameters from advanced parameters
#'   that typically do not require modification.
#' @param k Number of nearest neighbors used to construct the network.
#' @param th_anisotropy Threshold for anisotropy. Points exceeding this value are likely to be
#'   classified as wood.
#' @param min_passage Minimum number of times a point must be part of a least cost path to be
#'   considered part of the tree skeleton. Default: 5.
#' @param space_res Spatial resolution (in meters) for pathfinding. The algorithm does not actually
#'   evaluate "every point in space" as decribed above, but instead considers points at this interval.
#'   Default: 40 cm.
#' @param z_factor Numeric. Scaling factor for the Z-axis. The point cloud is compressed vertically to
#'   reduce the effect of large vertical gaps in the canopy while preserving horizontal distances.
#'   Default is 0.8.
#' @md
#' @export
segment_foliage = function(las, dtm, res = 0.05, max_gap = 1, th_anisotropy = 0.75, ...,  min_passage = 5, k = 5, space_res = 0.4, z_factor = 0.8)
{
  #res = 0.05; min_passage = 5; max_gap = 1;th_anisotropy = 0.75;k = 10; space_res = 0.4; z_factor = 0.8

  # Other parameters are hard coded
  connected_components_res = 0.05
  connected_components_min = 1000
  wood_assignation_k = 50
  wood_assignation_dist = 0.05
  wood_extra_reasignation_k = 10
  wood_extra_reasignation_dist = 0.03

  # The point cloud must have hag and anisotropy computed
  attributes <- names(las)
  stopifnot("anisotropy" %in% attributes)
  stopifnot("hag" %in% attributes)

  if (!"pointID" %in% names(las)) las@data$pointID = 1:lidR::npoints(las)

  . <- treeID <- X <- Y <- Z <-  hag <- hag_max <- hag_min <- anisotropy <- pointID <- wood <- decimated <- NULL

  ti <- tic()
  cat("Point cloud decimation... (1/11)\n") ; t0 = tic()

  # We will apply path finding to a decimated point cloud
  ans <- decimate_translate(las, res, z_factor)
  dec <- ans$dec
  dec@data <- dec@data[, .(X,Y,Z, anisotropy, pointID)]
  num_points <- lidR::npoints(dec)

  # Global translation to origin for computation stability
  # (translation already applied to dec)
  x_translation <- ans$x_translation
  y_translation <- ans$y_translation

  # The target points spread on the volume. These are end points
  # for the pathfinder
  target <- barycentric_predecimation(dec, space_res)
  target <- lidR::filter_poi(target, decimated == TRUE)
  target@data <- target@data[, .(X,Y,Z, anisotropy, pointID)]
  num_target <- lidR::npoints(target)

  # A layer of ground points used as a connector to the master seed
  gnd   <- seed_from_dtm(dtm)
  gnd$X <- gnd$X - x_translation
  gnd$Y <- gnd$Y - y_translation
  gnd$Z <- gnd$Z * z_factor
  lidR::quantize(gnd[["X"]], 0.01, las@header$`X offset`)
  lidR::quantize(gnd[["Y"]], 0.01, las@header$`Y offset`)
  lidR::quantize(gnd[["Z"]], 0.01, las@header$`Z offset`)
  gnd$anisotropy <- 1
  gnd$pointID    <- 0
  header  <- rlas::header_create(gnd)
  gnd     <- suppressWarnings(lidR::LAS(gnd, header))
  num_gnd <- lidR::npoints(gnd)

  # The master seed
  master_seed <- make_master_seed(gnd)

  # Plot for debugging
  if (FALSE)
  {
    x = plot(dec)
    plot(target, add = x, pal = "red", size = 4)
    plot(gnd, add = x, pal = "green", size = 6)
    plot(master_seed, add = x, pal = "white", size = 8)
  }

  toc(t0)
  cat("Building point cloud connectivity... (2/11)\n") ; t0 = tic()

  # Each point is connected to its knn. The connection is bidirectional. The cost of the connection
  # is based on the euclidean distance with some variation in order to specifically follow the wood
  # and not foliage. The cost is not the same in both directions

  point_network <- compute_point_network(dec, k = k, max_gap = max_gap)
  points_ids <- 1:num_points
  gc()

  # The cost is weighted by the anisotropy
  #A1 = point_coordinates$anisotropy[point_network$from]
  #A2 = point_coordinates$anisotropy[point_network$to]
  #W = 1-((A1+A2)/2)
  #W <- 0.1 + W * (0.9 - 0.1)
  #point_network$cost = point_network$cost * W
  #free(W, A1, A2)

  toc(t0)
  cat("Building target connectivity... (3/11)\n") ; t0 = tic()

  # Point cloud is connected to targets to reach. The connection is undirectional. It is possible to
  # move from the point cloud to the targets not the opposite. The cost of the connection is
  # null because the target are subsampled from the point cloud so distances are 0

  target_network <- compute_network(dec, target, k = 1)
  from <- target_network$from
  target_network$from <- target_network$to # switch direction
  target_network$to <- from
  target_network$to <- target_network$to + num_points
  target_ids <- (num_points+1):(num_points+1+num_target)

  toc(t0)
  cat("Building ground connectivity... (4/11)\n") ; t0 = tic()

  # The ground points are connected to their knn points in the scene. The connection is undirectional.
  # it is possible to move from the ground to the scene not the opposite. The cost of the connection
  # is the euclidean distance.

  ground_network <- compute_network(dec, gnd, k = k*10)
  ground_network$from <- ground_network$from + num_points + num_target
  ground_ids <- (num_points+num_target+1):(num_points+num_target+1+num_gnd)

  toc(t0)
  cat("Building master seed connectivity... (5/11)\n") ; t0 = tic()

  # The MASTER seed is connected to all the ground points in the scene. The connection cost is constant
  # and virtually 0. It allows to start from a single point to reach all the targets

  master_seed_network <- compute_network(gnd, master_seed, k = num_gnd)
  master_seed_network$from <- master_seed_network$from + max(ground_ids)
  master_seed_network$to <- master_seed_network$to + min(ground_ids)
  master_seed_network$cost <- 0.001
  master_seed_id <- master_seed_network[1,1]

  toc(t0)
  cat("Constructing the graph... (6/11)\n") ; t0 = tic()

  combined_network <- rbind(point_network, target_network, ground_network, master_seed_network)
  free(point_network, target_network, ground_network, master_seed_network)

  graph <- build_graph(combined_network)
  cache <- compute_distances(graph, master_seed_id)

  toc(t0)
  cat("Pathfinder... (7/11)\n") ; t0 = tic()

  # We maintain a counter for each point to count how many times the path find moved
  # by this point
  dec@data$count <- 0

  # We are moving from THE seed (repeated to match target size) to the targets
  from <- rep(master_seed_id, length(target_ids))
  to   <- target_ids

  # For loop by chunk to reduce memory usage and have an estimated progression
  chunk_size <- 50000
  chunks     <- split(to, ceiling(seq_along(to) / chunk_size))
  pb <- utils::txtProgressBar(min = 0, max = length(chunks), style = 3, width = 50)

  for (i in seq_along(chunks))
  {
    current_to <- chunks[[i]]

    path <- findPaths(graph, cache, from[seq_along(current_to)], current_to)
    path <- path$paths
    path <- lapply(path, function(x) x[-1])
    path <- do.call(c, path)

    count <- table(path)
    id    <- as.numeric(names(count))

    rm    <- id > num_points
    id    <- id[!rm]
    count <- count[!rm]

    dec@data$count[id] <- dec@data$count[id] + count

    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  free(combined_network)

  # We now know, for each point, how many times the pathfinder moved by this points
  dec@data$count[dec@data$count > 0] <- log(dec@data$count[dec@data$count > 0])

  free(path, count, id, rm, from, to)

  toc(t0)
  cat("Assigning wood to small structure... (8/11)\n") ; t0 = tic()

  skeleton <- lidR::filter_poi(dec, count > log(min_passage))

  las@data$passage <- 0
  las@data$passage[skeleton$pointID] <- exp(skeleton$count)
  las <- lidR::add_lasattribute_manual(las, name = "passage", desc = "passage points", type = "int")

  # The decimated points
  skeleton$X <- skeleton$X + x_translation
  skeleton$Y <- skeleton$Y + y_translation
  skeleton$Z <- skeleton$Z / z_factor

  skeleton_neighbors  <- lidR::knnx(las, skeleton, k = wood_assignation_k)
  rm <- skeleton_neighbors$nn.dist > wood_assignation_dist
  id <- skeleton_neighbors$nn.index[!rm]

  las@data$wood = FALSE
  las@data$wood[id] = TRUE

  free(skeleton_neighbors, rm, id)

  toc(t0)
  cat("Filter anisotropy... (Step 9/11)\n") ; t0 = tic()

  # Remove foliage based on anisotropy only

  nofoliage = lidR::filter_poi(las, anisotropy > th_anisotropy | wood == TRUE)

  toc(t0)
  cat("Connected component cleaning... (10/11)\n") ; t0 = tic()

  # Looking only at wood points (high anisotropy), perform a connected component
  # scan, remove patches with not enough points. They are reasigned as foliage

  nofoliage$Z <- nofoliage$Z * 0.8
  nofoliage   <- lidR::connected_components(nofoliage, connected_components_res, connected_components_min)
  nofoliage   <- nofoliage[nofoliage$clusterID != 0 | nofoliage$wood == TRUE]
  nofoliage$Z <- nofoliage$Z / 0.8

  toc(t0)
  cat("Extra wood reasignation... (11/11)\n") ; t0 = tic()

  # We look at the neighboring points of the wood.  Points close to the wood
  # are wood points too. This assigns extra wood point is the branches and remove
  # some false negatives

  wood_neighbors <- lidR::knnx(las, nofoliage, k = wood_extra_reasignation_k)
  rm <- wood_neighbors$nn.dist > wood_extra_reasignation_dist
  id <- wood_neighbors$nn.index[!rm]

  las@data$foliage <- TRUE
  las@data$foliage[id] <- FALSE
  las@data$foliage <- as.integer(las$foliage)
  las <- lidR::add_lasattribute_manual(las, name = "foliage", desc = "foliage 1 wood 0", type = "char")

  free(nofoliage, rm, id, wood_neighbors)

  toc(t0)
  toc(ti, space = "")
  gc()

  return(las)
}

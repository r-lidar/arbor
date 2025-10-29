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
#' @param params list See \link{parameters}.
#' @md
#' @export
segment_foliage = function(las, dtm, params = default_parameters)
{
  #res = 0.05; min_passage = 5; max_gap = 1;th_anisotropy = 0.75;k = 10; space_res = 0.4; z_factor = 0.8

  z_factor <- params$semantic$z_scale

  # The point cloud must have hag and anisotropy computed
  attributes <- names(las)
  stopifnot("anisotropy" %in% attributes)
  stopifnot("hag" %in% attributes)

  if (!"pointID" %in% names(las)) las@data$pointID = 1:lidR::npoints(las)

  . <- treeID <- X <- Y <- Z <-  hag <- hag_max <- hag_min <- anisotropy <- pointID <- wood <- decimated <- NULL

  ti <- tic()
  cat("Point cloud decimation... (1/12)\n") ; t0 = tic()

  # We will apply path finding to a decimated point cloud
  ans        <- decimate_translate(las, params)
  dec        <- ans$dec
  dec@data   <- dec@data[, .(X,Y,Z, anisotropy, pointID)]
  dec$Z      <- dec$Z * z_factor

  # Global translation to origin for computation stability
  # (translation already applied to dec)
  x_translation <- ans$x_translation
  y_translation <- ans$y_translation

  # The target points spread on the volume. These are end points
  # for the pathfinder
  target      <- barycentric_decimation(dec, params$path_finder$res)
  target@data <- target@data[, .(X,Y,Z, anisotropy, pointID)]

  # A layer of ground points used as a connector to the master seed
  gnd <- make_ground_points(dtm, x_translation, y_translation, z_factor, params$semantic$ground_res)

  # The master seed
  master_seed <- make_master_seed(gnd)

  # Plot for debugging
  if (FALSE)
  {
    x <- plot(dec)
    plot(target, add = x, pal = "red", size = 4)
    plot(gnd, add = x, pal = "green", size = 6)
    plot(master_seed, add = x, pal = "white", size = 8)
  }

  toc(t0)
  cat("Building point cloud connectivity... (2/12)\n") ; t0 = tic()

  # Each point is connected to its knn. The connection is bidirectional. The cost of the connection
  # is based on the euclidean distance with some variation in order to specifically follow the wood
  # and not foliage. The cost is not the same in both directions
  k <- params$path_finder$k_neighborhood_connectivity
  max_gap <- params$path_finder$max_gap

  graph <- build_semantic_graph(dec@data, target@data, gnd@data, master_seed@data, k, max_gap)

  # The cost is weighted by the anisotropy
  #A1 <- dec$anisotropy[point_network$from]
  #A2 <- dec$anisotropy[point_network$to]
  #W  <- 1-(A1+A2)/2
  #point_network$cost = point_network$cost * W
  #free(W, A1, A2)

  toc(t0)
  cat("Pathfinder... (7/12)\n") ; t0 = tic()

  num_points <- lidR::npoints(dec)
  num_target <- lidR::npoints(target)
  num_gnd    <- lidR::npoints(gnd)
  target_ids <- 1:num_target + num_points - 1
  master_id  <- num_points + num_target + num_gnd
  dec@data$passage <- accumulate_passages(graph, master_id, target_ids, num_points)

  if (FALSE) plot_passage(dec)

  free(graph)

  las@data$passage <- 0
  las@data$passage[dec$pointID] <- dec$passage
  las <- lidR::add_lasattribute_manual(las, name = "passage", desc = "passage points", type = "int")
  gc()

  toc(t0)
  cat("Assigning wood to small structure... (8/12)\n") ; t0 = tic()

  # The decimated points
  min_passage <- params$path_finder$min_passage
  wood_assignation_k <- params$semantic$wood_assignation_k
  wood_assignation_dist <- params$semantic$wood_assignation_dist

  passage    <- NULL
  skeleton   <- lidR::filter_poi(dec, passage > min_passage)
  skeleton$X <- skeleton$X + x_translation
  skeleton$Y <- skeleton$Y + y_translation
  skeleton$Z <- skeleton$Z / z_factor # Rescale true Z to rematch original point cloud

  skeleton_neighbors  <- lidR::knnx(las, skeleton, k = wood_assignation_k)
  rm <- skeleton_neighbors$nn.dist > wood_assignation_dist
  id <- skeleton_neighbors$nn.index[!rm]

  # 100% sure those ones are wood
  path_finder_based_wood <- rep(FALSE, lidR::npoints(las))
  path_finder_based_wood[id] <- TRUE

  las@data$wood <- FALSE
  las@data$wood[id] <- TRUE

  # Plot for debugging
  if (FALSE) plot(las, color = "wood", pal = rev(foliage.colors))

  free(skeleton_neighbors, rm, id)

  toc(t0)
  cat("Filter high anisotropy... (9/12)\n") ; t0 = tic()

  # Remove foliage based on anisotropy only
  th_high_  <- params$semantic$high_anisotropy_threshold
  connected_components_res <- params$semantic$connected_components_res
  connected_components_min <- params$semantic$connected_components_min

  nofoliage <- lidR::filter_poi(las, anisotropy > th_high_ | wood == TRUE)

  # Plot for debuging
  if (FALSE)
  {
    foliage <- lidR::filter_poi(las, anisotropy < 0.75)
    foliage <- compute_anisotropy(foliage, k = 50)
    plot(foliage, color = "anisotropy", legend = T, breaks = "quantile")

    plot(foliage)
    plot(nofoliage, pal = foliage.colors[1])
    plot(nofoliage, color = "anisotropy", legend = T, breaks = "quantile")
    plot(foliage, pal = foliage.colors[2])
  }

  # Looking only at wood points (high anisotropy), perform a connected component
  # scan, remove patches with not enough points. They are reasigned as foliage

  nofoliage$Z <- nofoliage$Z * z_factor
  nofoliage   <- lidR::connected_components(nofoliage, connected_components_res, connected_components_min)
  nofoliage$clusterID[nofoliage$clusterID == 0] = NA_integer_
  if (FALSE) plot(nofoliage, color = "clusterID") # Plot for debuging
  nofoliage   <- nofoliage[!is.na(nofoliage$clusterID)]
  nofoliage   <- nofoliage[nofoliage$clusterID != 0 | nofoliage$wood == TRUE]
  nofoliage$Z <- nofoliage$Z / z_factor


  if (FALSE) plot(nofoliage) # Plot for debugging

  # 100% sure those ones are wood
  high_anisotropy_wood <- rep(FALSE, lidR::npoints(las))
  high_anisotropy_wood[nofoliage$pointID] <- TRUE

  toc(t0)
  cat("Filter medium anisotropy... (Step 10/12)\n") ; t0 = tic()

  th_medium_ <- params$semantic$medium_anisotropy_thresold
  sor_k      <- params$semantic$medium_anisotropy_sor_k
  sor_m      <- params$semantic$medium_anisotropy_sor_m
  connected_components_res <- params$semantic$connected_components_res
  connected_components_min <- params$semantic$connected_components_min

  nofoliage <- lidR::filter_poi(las, (anisotropy > th_medium_ & anisotropy < th_high_)  | wood == TRUE)
  nofoliage <- lidR::classify_noise(nofoliage, lidR::sor(sor_k, sor_m))
  if (FALSE) plot(nofoliage, color = "Classification") # Plot for debuging
  nofoliage <- lidR::remove_noise(nofoliage)
  nofoliage$Z <- nofoliage$Z * z_factor
  nofoliage   <- lidR::connected_components(nofoliage, connected_components_res, connected_components_min)
  nofoliage$clusterID[nofoliage$clusterID == 0] <- NA_integer_
  if (FALSE) plot(nofoliage, color = "clusterID") # Plot for debuging
  nofoliage   <- nofoliage[!is.na(nofoliage$clusterID)]
  nofoliage$Z <- nofoliage$Z / z_factor

  # Those ones are likely wood
  medium_anistropy_wood <- rep(FALSE, lidR::npoints(las))
  medium_anistropy_wood[nofoliage$pointID] <- TRUE

  las@data$wood <- path_finder_based_wood | medium_anistropy_wood | high_anisotropy_wood

  if (FALSE) plot(las, color = "wood", pal = rev(foliage.colors)) # Plot for debuging

  toc(t0)
  cat("Extra wood reasignation... (11/12)\n") ; t0 = tic()

  # We look at the neighboring points of the wood.  Points close to the wood
  # are wood points too. This assigns extra wood point is the branches and remove
  # some false negatives
  wood_extra_reasignation_k    <- params$semantic$wood_extra_reasignation_k
  wood_extra_reasignation_dist <- params$semantic$wood_extra_reasignation_dist

  nofoliage <- lidR::filter_poi(las, wood == TRUE)

  wood_neighbors <- lidR::knnx(las, nofoliage, k = wood_extra_reasignation_k)
  rm <- wood_neighbors$nn.dist > wood_extra_reasignation_dist
  id <- wood_neighbors$nn.index[!rm]

  las@data$foliage <- 1
  las@data$foliage[id] <- 0
  las <- lidR::add_lasattribute_manual(las, name = "foliage", desc = "foliage: 1 or 2 wood: 0", type = "char")

  if (FALSE) plot(las, color = "foliage", pal = foliage.colors) # Plot for debuging

  free(nofoliage, rm, id, wood_neighbors)

  toc(t0)
  cat("Extra foliage reasignation... (12/12)\n") ; t0 = tic()

  foliage <- lidR::filter_poi(las, foliage == 1)
  if (FALSE) plot(foliage, color = "anisotropy", legend = T, breaks = "quantile")
  high_ani <- lidR::filter_poi(foliage, anisotropy > th_high_)
  las@data$foliage[high_ani$pointID] <- 2

  free(high_ani)

  toc(t0)
  toc(ti, space = "")
  gc()

  return(las)
}

make_ground_points = function(dtm, xoffset, yoffset, zscale, res)
{
  gnd   <- seed_from_dtm(dtm, res = res)
  gnd$X <- gnd$X - xoffset
  gnd$Y <- gnd$Y - yoffset
  gnd$Z <- gnd$Z * zscale
  lidR::quantize(gnd[["X"]], 0.01, las@header[["X offset"]])
  lidR::quantize(gnd[["Y"]], 0.01, las@header[["Y offset"]])
  lidR::quantize(gnd[["Z"]], 0.01, las@header[["Z offset"]])
  header  <- rlas::header_create(gnd)
  gnd     <- suppressWarnings(lidR::LAS(gnd, header))
  gnd
}

seed_from_dtm = function(dtm, res)
{
  seeds = terra::rast(terra::ext(dtm), res = res)
  seeds = terra::resample(dtm, seeds)
  seeds = as.data.frame(seeds, xy = T)
  seeds = data.table::as.data.table(seeds)
  names(seeds) = c("X", "Y", "Z")
  seeds
}

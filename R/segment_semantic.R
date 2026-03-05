#' Wood-Foliage Semantic Segmentation
#'
#' This function performs wood-foliage semantic segmentation by constructing a network of points
#' using k-nearest neighbors (KNN). A pathfinding algorithm navigates this network to identify the
#' tree skeleton by finding the least-cost path from every point in space to the ground. Each point
#' is classified as wood or foliage based on two main criteria: (1) Proximity to the detected skeleton and
#' (2) wood likelihood exceeding a given threshold. The function \link{wood_likelihood} must be applied
#' first. The point cloud must have an attribute 'hag' (see \link[lidR:height_above_ground]{height_above_ground})
#' Finally, a connected component step removes small clusters incorrectly classified as wood.
#'
#' @param las A LAS object from lidR.
#' @param dtm A SpatRaster object. Digital Terrain Model
#' @param params list See \link{parameters}.
#' @md
#' @rdname segment_semantic
#' @name segment_semantic
#' @export
segment_semantic = function(las, dtm, params = default_arbor_parameters)
{
  params <- evaluate_penalty(params)
  las@data$passage <- integer(lidR::npoints(las))
  las@data$foliage <- integer(lidR::npoints(las))
  las <- lidR::add_lasattribute_manual(las, name = "passage", desc = "passage points", type = "int")
  las <- lidR::add_lasattribute_manual(las, name = "foliage", desc = "foliage: 1 or 2 wood: 0", type = "char")
  gnd <- make_ground_points(dtm, params$semantic$ground_res, las@header)
  segment_semantic_cpp(las@data, gnd@data, params)
  return(las)
}

make_ground_points = function(dtm, res, header)
{
  gnd   <- seed_from_dtm(dtm, res = res)
  lidR::quantize(gnd[["X"]], 0.01, header[["X offset"]])
  lidR::quantize(gnd[["Y"]], 0.01, header[["Y offset"]])
  lidR::quantize(gnd[["Z"]], 0.01, header[["Z offset"]])
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

# # Old R around C++
# segment_semantic_r_cpp = function(las, dtm, params = default_arbor_parameters)
# {
#   stopifnot("pwood" %in% names(las))
#
#   logger("Segmentic segmentation start")
#
#   # Accumulate passages
#   params   <- evaluate_penalty(params)
#   gnd      <- make_ground_points(dtm, params$semantic$ground_res, las@header)
#   passages <- accumulate_passages_cpp(las@data, gnd@data, params)
#   las      <- lidR::add_lasattribute_manual(las, passages, name = "passage", desc = "passage points", type = "int")
#
#   # Wood from pathfinder
#   path_finder_based_wood = assign_wood_from_passage_cpp(las@data, params)
#   las@data$foliage = as.integer(!path_finder_based_wood)
#
#   # Wood from likelihood
#   high_likelihood_based_wood = assign_wood_from_high_likelihood_cpp(las@data, params)
#   medium_likelihood_based_wood = assign_wood_from_medium_likelihood_cpp(las@data, params)
#
#   is_wood = path_finder_based_wood | high_likelihood_based_wood | medium_likelihood_based_wood
#   las@data$foliage = as.integer(!is_wood)
#
#   # Wood from dilatation
#   is_wood = assign_wood_from_wood_dilatation_cpp(las@data, params)
#   las@data$foliage = as.integer(!is_wood)
#
#   logger("Extra class 2 foliage re-assignation... (9/10)")
#   las@data[foliage == 1 & pwood > params$semantic$high_pwood_threshold, foliage := 2L]
#
#   las <- lidR::add_lasattribute_manual(las, name = "foliage", desc = "foliage: 1 or 2 wood: 0", type = "char")
#
#   logger("Semantic segmentation completed")
#   gc()
#
#   return(las)
# }

# Old R code for debugging and rendering
# segment_semantic_r = function(las, dtm, params = default_arbor_parameters)
# {
#   attributes <- names(las)
#   stopifnot("pwood" %in% attributes)
#   stopifnot("hag" %in% attributes)
#   . <- treeID <- X <- Y <- Z <-  hag <- hag_max <- hag_min <- pwood <- pointID <- wood <- decimated <- NULL
#
#   if (!"pointID" %in% names(las)) las@data$pointID = 1:lidR::npoints(las)
#
#   ti <- tic() ; cat("Point cloud decimation... (1/8)\n") ; t0 = tic()
#
#   core   <- hybrid_homogeneization(las, params$path_finder$decimation)
#   target <- barycentric_decimation(las, params$path_finder$space_res)
#   gnd    <- make_ground_points(dtm, params$semantic$ground_res, las@header)
#
#   # Plot for debugging
#   if (FALSE)
#   {
#     x <- plot(core)
#     plot(target, add = x, pal = "red", size = 4)
#     plot(gnd, add = x, pal = "green", size = 6)
#   }
#
#   toc(t0) ; cat("Building point cloud connectivity... (2/8)\n") ; t0 = tic()
#
#   params <- evaluate_penalty(params)
#   graph  <- build_semantic_graph(core@data, target@data, gnd@data, params)
#
#   # The cost is weighted by the pwood
#   #A1 <- core$pwood[point_network$from]
#   #A2 <- core$pwood[point_network$to]
#   #W  <- 1-(A1+A2)/2
#   #point_network$cost = point_network$cost * W
#   #free(W, A1, A2)
#
#   toc(t0); cat("Pathfinder... (3/8)\n") ; t0 = tic()
#
#   num_points <- lidR::npoints(core)
#   num_target <- lidR::npoints(target)
#   num_gnd    <- lidR::npoints(gnd)
#   target_ids <- 1:num_target + num_points - 1
#   ground_ids <- 1:num_gnd + num_target + num_points - 1
#   master_id  <- num_points + num_target + num_gnd
#
#   core@data$passage <- accumulate_passages_old(graph, master_id, target_ids, num_points)
#
#   if (FALSE)
#   {
#     x = plot(core, pal = "gray")
#     plot_passage(core, add = x, size = 3)
#     plot(gnd, add = x, pal = "darkgreen", size = 3)
#   }
#
#   free(graph)
#
#   las@data$passage <- 0
#   las@data$passage[core$pointID] <- core$passage
#   las <- lidR::add_lasattribute_manual(las, name = "passage", desc = "passage points", type = "int")
#
#   toc(t0) ; cat("Assigning wood to small structure... (4/8)\n") ; t0 = tic()
#
#   min_passage           <- params$semantic$min_passage
#   wood_assignation_k    <- 10
#   wood_assignation_dist <- params$semantic$wood_assignation_dist
#
#   passage    <- NULL
#   skeleton   <- lidR::filter_poi(core, passage > min_passage)
#
#   skeleton_neighbors  <- lidR::knnx(las, skeleton, k = wood_assignation_k)
#   rm <- skeleton_neighbors$nn.dist > wood_assignation_dist
#   id <- skeleton_neighbors$nn.index[!rm]
#
#   # 100% sure those ones are wood because on the pathfinder
#   path_finder_based_wood <- rep(FALSE, lidR::npoints(las))
#   path_finder_based_wood[id] <- TRUE
#
#   las@data$wood <- FALSE
#   las@data$wood[id] <- TRUE
#
#   if (FALSE) plot(las, color = "wood", pal = rev(foliage.colors)) # Plot for debugging
#
#   free(skeleton_neighbors, rm, id)
#
#   toc(t0) ; cat("Filter high likeliwood... (5/8)\n") ; t0 = tic()
#
#   z_factor <- 1
#
#   # Remove foliage based on high likelihood only. High anistropy = wood
#
#   th_high_                 <- params$semantic$high_pwood_threshold
#   connected_components_res <- params$semantic$connected_components_res
#   connected_components_min <- params$semantic$connected_components_min
#
#   nofoliage <- lidR::filter_poi(las, pwood > th_high_ | wood == TRUE)
#
#   # Plot for debuging
#   if (FALSE)
#   {
#     plot(nofoliage, pal = foliage.colors[1])
#     plot_likelihood(nofoliage)
#   }
#
#   # Looking only at wood points (high likelihood), perform a connected component
#   # scan, remove patches with not enough points. They are reasigned as foliage
#
#   nofoliage$Z <- nofoliage$Z * z_factor
#   nofoliage   <- connected_components(nofoliage, connected_components_res, connected_components_min, connectivity = 26)
#   nofoliage$clusterID[nofoliage$clusterID == 0] = NA_integer_
#   if (FALSE) plot(nofoliage, color = "clusterID") # Plot for debuging
#   nofoliage   <- nofoliage[!is.na(nofoliage$clusterID)]
#   nofoliage   <- nofoliage[nofoliage$clusterID != 0 | nofoliage$wood == TRUE]
#   nofoliage$Z <- nofoliage$Z / z_factor
#
#   if (FALSE) plot(nofoliage) # Plot for debugging
#
#   # 100% sure those ones are wood
#   high_pwood_wood <- rep(FALSE, lidR::npoints(las))
#   high_pwood_wood[nofoliage$pointID] <- TRUE
#
#   toc(t0) ; cat("Filter medium likelihood (Step 6/8)\n") ; t0 = tic()
#
#   th_medium_               <- params$semantic$medium_pwood_thresold
#   sor_k                    <- params$semantic$medium_pwood_sor_k
#   sor_m                    <- params$semantic$medium_pwood_sor_m
#   connected_components_res <- params$semantic$connected_components_res
#   connected_components_min <- params$semantic$connected_components_min
#
#   nofoliage <- lidR::filter_poi(las, (pwood > th_medium_ & pwood < th_high_)  | wood == TRUE)
#   nofoliage <- sor(nofoliage, sor_k, sor_m)
#   if (FALSE) plot(nofoliage, color = "Classification") # Plot for debuging
#   nofoliage <- lidR::remove_noise(nofoliage)
#   nofoliage$Z <- nofoliage$Z * z_factor
#   nofoliage   <- connected_components(nofoliage, connected_components_res, connected_components_min, connectivity = 26)
#   nofoliage$clusterID[nofoliage$clusterID == 0] <- NA_integer_
#   if (FALSE) plot(nofoliage, color = "clusterID") # Plot for debuging
#   nofoliage   <- nofoliage[!is.na(nofoliage$clusterID)]
#   nofoliage$Z <- nofoliage$Z / z_factor
#
#   # Those ones are likely wood
#   medium_anistropy_wood <- rep(FALSE, lidR::npoints(las))
#   medium_anistropy_wood[nofoliage$pointID] <- TRUE
#
#   las@data$wood <- path_finder_based_wood | medium_anistropy_wood | high_pwood_wood
#
#   if (FALSE) plot(las, color = "wood", pal = rev(foliage.colors)) # Plot for debuging
#
#   toc(t0)
#   cat("Extra wood reasignation... (7/8)\n") ; t0 = tic()
#
#   # We look at the neighboring points of the wood.  Points close to the wood
#   # are wood points too. This assigns extra wood point is the branches and remove
#   # some false negatives
#   wood_extra_reasignation_k    <- params$semantic$wood_extra_reasignation_k
#   wood_extra_reasignation_dist <- params$semantic$wood_extra_reasignation_dist
#
#   nofoliage <- lidR::filter_poi(las, wood == TRUE)
#
#   wood_neighbors <- lidR::knnx(las, nofoliage, k = wood_extra_reasignation_k)
#   rm <- wood_neighbors$nn.dist > wood_extra_reasignation_dist
#   id <- wood_neighbors$nn.index[!rm]
#
#   las@data$foliage <- 1
#   las@data$foliage[id] <- 0
#   las <- lidR::add_lasattribute_manual(las, name = "foliage", desc = "foliage: 1 or 2 wood: 0", type = "char")
#
#   if (FALSE) plot(las, color = "foliage", pal = foliage.colors) # Plot for debuging
#
#   free(nofoliage, rm, id, wood_neighbors)
#
#   toc(t0) ; cat("Extra foliage reasignation... (8/8)\n") ; t0 = tic()
#
#   foliage <- lidR::filter_poi(las, foliage == 1)
#   if (FALSE) plot_likelihood(foliage)
#   high_ani <- lidR::filter_poi(foliage, pwood > th_high_)
#   las@data$foliage[high_ani$pointID] <- 2
#
#   free(high_ani)
#
#   toc(t0)
#   toc(ti, space = "")
#   gc()
#
#   return(las)
# }

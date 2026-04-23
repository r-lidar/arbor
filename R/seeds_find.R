#' Find seeds to perform instance segmentation
#'
#' In order to perform instance segmentation with  \link{segment_instance} we need
#' some seeds with reference treeIDs. Thus function finds the seeds
#'
#' @param las A LAS object from lidR.
#' @param params list See \link{parameters}.
#' @export
find_seeds <- function(las, params)
{
  #params = evaluate_penalty(params)
  cloud = las@data ; cloud[["treeID"]] = -1L # Allocate memory because C++ assumes memory allocated
  seeds = find_seeds_cpp(cloud, params)
  seeds = suppressWarnings(lidR::LAS(seeds, lidR::header(las)))
  return(seeds)
}

# find_seeds_r <- function(las, params)
# {
#   # The point cloud is supposed to have passage hag and foliage
#   attributes <- names(las)
#   stopifnot("hag" %in% attributes)
#   stopifnot("passage" %in% attributes)
#   stopifnot("foliage" %in% attributes)
#
#   logger("Seed detection start")
#
#   logger("Extracting passages");
#   long_passages  <- extract_passages(las, "long", params)
#   short_passages <- extract_passages(las, "short", params)
#   if (FALSE) { x = plot(long_passages, pal = "green") ; plot(short_passages, add = x, pal = "red") }
#
#   logger("Slice wood")
#   wood <- slice_wood(las, heights, thick)
#   if (FALSE) { x = plot(wood) ; plot(long_passages, add = x, pal = "green") }
#
#
#   logger("Circle detection")
#   circles = detect_tree_circles(wood) # r_detect_tree_circles(wood)
#   if (nrow(circles) == 0) stop("No circle detected in wood slices")
#
#   if (FALSE) {
#     pal = pastel.colors(max(circles$id))
#     col = pal[circles$id]
#     x <- plot(wood)
#     plot(long_passages, pal = "green", add = x)
#     for (i in 1:nrow(circles)) add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], circles$Z[i], col = col[i])
#   }
#
#
#   logger("Generate tree cages")
#   cage = generate_cage(circles, params, lidR::header(las))
#   if (FALSE) { plot(cage) }
#
#
#   logger("Safe zone exclusion")
#   wood = safe_zone(wood, circles)
#   if (FALSE) { x = plot(wood) ; plot(long_passages, add = x, pal = "green") }
#
#
#   logger("Find primary seeds")
#   primary_seeds = find_primary_seeds(wood, long_passages, cage) # We have seed for the big trees (long passage)
#   if (FALSE) { x <- plot(primary_seeds, color = "treeID") ; plot(short_passages, add = x, pal = "red")}
#
#
#   logger("Pathfinder for minor tree seeds")   # Some short passage could be from big tree. Path finder to attach them.
#   short_passages = merge_short_passages(short_passages, primary_seeds)
#   short_passages_withid = short_passages[!is.na(short_passages$treeID)]
#   short_passages_withid$foliage = NULL
#   primary_seeds = rbind(primary_seeds, short_passages_withid)
#   if (FALSE) { plot(primary_seeds, color = "treeID") }
#
#
#   logger("Find secondary seeds")
#   max_id = max(seeds$treeID, na.rm = TRUE)
#   secondary_seeds = find_secondary_seeds(short_passages, max_id)
#   if (FALSE) { x = plot(primary_seeds, color = "treeID") ; plot(secondary_seeds, add = x, color = "treeID") }
#
#   logger("Merge primary and secondary seeds")
#   seeds = rbind(primary_seeds, secondary_seeds)
#   seeds = lidR::filter_poi(seeds, passage > 0, hag < 1)
#   if (FALSE) { plot(seeds, color = "treeID") }
#
#
#   logger("Seed detection completed")
#
#   return(seeds)
# }
#
# extract_passages = function(las, type, params)
# {
#   X <- Y <- Z <- passage <- hag <- NULL
#
#   # The heights at which we slice
#   # ( we extract some slices of wood (thickness 3cm))
#   thick   <- params$seed$slice_thickness
#   heights <- c(min(las$hag)+thick, params$seed$slice_at)
#
#   # Extract the passages of big trees and small features
#   th                   <- max(heights) + 0.2
#   min_passage          <- params$seed$min_passage
#
#   if (type == "long")
#   {
#     long_passages        <- lidR::filter_poi(las, passage > min_passage, hag < th)
#     long_passages@data   <- long_passages@data[, .(X,Y,Z, passage, hag)]
#     long_passages@data   <- densify_passage(long_passages@data)
#     return(long_passages)
#   }
#
#   if (type == "short")
#   {
#     short_passages       <- lidR::filter_poi(las, passage > 1, hag < min(hag) + 0.5, passage < 10)
#     short_passages@data  <- short_passages@data[, .(X,Y,Z, passage, hag)]
#     short_passages@data  <- densify_passage(short_passages@data)
#     return(short_passages)
#   }
#
#   stop("Internal error: invalid type")
# }
#
# slice_wood = function(las, heights, thinkness = 0.02)
# {
#   hag <- X <- Y <- Z <- . <- passage <- NULL
#
#   logger("Slicing the point cloud")
#
#   slice_filter <- Reduce(`|`, lapply(heights, function(s)
#   {
#     (las$hag > (s-thinkness/2) & las$hag < (s + thinkness/2))
#   }))
#
#   wood   <- lidR::filter_poi(las, slice_filter)
#   wood   <- lidR::filter_poi(wood, foliage == 0)
#   wood@data <- wood@data[, .(X,Y,Z, passage, hag)]
#
#   return(wood)
# }
#
# safe_zone = function(wood, circles)
# {
#   # For each circle we exclude wood in a safe zone beyond the circles.
#   # This allow to clean false positive around important trees and prevent
#   # dummy connection caused by noise in the next connected component step
#   px <- wood$X
#   py <- wood$Y
#   pz <- wood$Z
#   rm <- rep(FALSE, lidR::npoints(wood))
#   safe_zone <- 0.2
#   for (i in 1:nrow(circles))
#   {
#     cx <- circles$X[i]
#     cy <- circles$Y[i]
#     cz <- circles$Z[i]
#     r  <- circles$R[i]
#     d  <- sqrt((px-cx)^2 + (py-cy)^2 +(pz-cz)^2)
#     rm[d > (r + 0.02) & d  < (r + safe_zone)] <- TRUE
#   }
#   wood <- wood[!rm]
#
#   if (FALSE)
#   {
#     plot(wood)
#     for (i in 1:nrow(circles))
#       add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], circles$Z[i])
#   }
#
#   if (FALSE)
#   {
#     x = plot(long_passages, pal = "green")
#     plot(short_passages, add = x, pal = "red")
#     plot(wood, add = x, pal = foliage.colors[1])
#     for (i in 1:nrow(circles))
#       add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], circles$Z[i])
#   }
#
#   return(wood)
# }
#
# densify_passage <- function(data, offset = 0.01)
# {
#   data_up <- data
#   data_up[["Z"]] <- data[["Z"]] + offset
#   data_up[["passage"]] = -1
#
#   data_down <- data
#   data_down[["Z"]] <- data[["Z"]] - offset
#   data_down[["passage"]] = -1
#
#   return(rbind(data, data_up, data_down))
# }
#
# find_primary_seeds = function(wood, long_passages, cage)
# {
#   # Bind the wood, the long passages and the cage and compute connected component and merge passage from the same trees
#   lidR::st_crs(long_passages) = lidR::st_crs(wood)
#   lidR::st_crs(cage) = lidR::st_crs(wood)
#   temp   <- suppressWarnings(rbind(wood, long_passages, cage))
#   res    <- round(params$path_finder$decimation*0.8, 2)
#   temp$Z <- temp$Z * 0.5
#   temp   <- connected_components(temp, res, 1, name = "treeID", connectivity = 26)
#   temp$Z <- temp$Z / 0.5
#
#   if (FALSE) plot(temp, color = "treeID")
#
#   # We have seed for the big trees (long passage)
#   long_passages_seeds = lidR::filter_poi(temp, passage > 0)
#
#   long_passages_seeds
# }
#
# find_secondary_seeds = function(short_passages, id_start)
# {
#   # Some short passage don't have IDs. Assign an ID
#   short_passages_noid = short_passages[is.na(short_passages$treeID)]
#   short_passages_noid = connected_components(short_passages_noid, 0.1, 1, "treeID")
#   short_passages_noid$treeID = short_passages_noid$treeID + id_start + 1
#   short_passages_noid@data$foliage = NULL
#   return(short_passages_noid)
# }
#
# merge_short_passages = function(short_passages, long_passages_seeds)
# {
#   short_passages@data$foliage = 0 # Everything is wood
#   p = default_arbor_parameters
#   p$path_finder$max_gap = 0.1
#   p$path_finder$k_neighborhood_connectivity = 10
#   p$path_finder$k_seed_connectivity = 2
#   p$path_finder$distance_power = 1
#   p$path_finder$angle_penalty = function(x) {return(rep(1, length(x))) }
#
#   sink(tempfile())
#   on.exit(suppressWarnings(sink()), add = TRUE)
#   short_passages = segment_instance(short_passages, long_passages_seeds, p)
#   sink()
#
#   return(short_passages)
# }
#
# add_circle3d <- function(x, center_x, center_y, radius, height, col = "red")
# {
#   theta <- seq(0, 2 * pi, length.out = 50)
#   xx <- center_x - x[1] + radius * cos(theta)
#   yy <- center_y - x[2] + radius * sin(theta)
#   zz <- rep(height, 50)
#
#   # Plot the circle in 3D
#   rgl::lines3d(xx, yy, zz, lwd = 5, col = col)
# }
#
# qpoints3d = function(x)
# {
#   dx = mean(x[[1]])
#   dy = mean(x[[2]])
#   x[[1]] = x[[1]] - dx
#   x[[2]] = x[[2]] - dy
#   rgl::points3d(x)
# }


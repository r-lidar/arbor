
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
#   p = arbor_parameters_default
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


detect_tree_circles = function(wood)
{
  detect_tree_circles_cpp(wood@data)
}

# r_detect_tree_circles = function(wood)
# {
#
#   logger("  Connected component")
#
#   # Connect the wood point into clusters
#   cl_wood <- connected_components(wood, 0.05, 10, connectivity = 26)
#   cl_wood <- lidR::filter_poi(cl_wood, clusterID != 0)
#
#   if (FALSE)
#   {
#     x = plot(cl_wood, color = "clusterID", pal = pastel.colors(200))
#     plot(long_passages, add = x, pal = "green")
#   }
#
#   logger("  Circle detection per component")
#
#   # For each cluster search for circles. If we have a nice circle we have a tree
#   is.valid.circle <- function(radius, angle_range, pinliner, pinside)
#   {
#     if (radius > 2) return(FALSE)
#     if (radius < 0.02)  return(FALSE)
#     if (radius  < 0.05)  return(angle_range > 180 & pinliner > 30)
#     if (pinside > 20)   return(FALSE)
#     if (radius < 0.10)  return(angle_range > 130 & pinliner > 60)
#     return(angle_range > 140 & pinliner > 40)
#   }
#   fit_circle_to_seed <- function(cl)
#   {
#     id = cl$clusterID[1]
#     if (nrow(cl) < 20) return(NULL)
#     cl <- as.matrix(cl[,1:3])
#     circle <- ransac_circle(cl, num_iterations = 400, inlier_threshold = 0.02)
#
#     valid  <- is.valid.circle(circle$radius, circle$covered_arc_degree, circle$percentage_inlier*100, circle$percentage_inside*100)
#
#     if (FALSE)
#     {
#       if (valid) col = "darkgreen" else col = "red"
#       plot(cl[,1], cl[,2], asp = 1, main = paste(i, "id =", id))
#       inl = circle$inliers
#       points(cl[inl,1], cl[inl,2], pch = 19)
#       symbols(circle$center_x, circle$center_y, circles = circle$radius, inches = FALSE, add = TRUE, fg = col)
#       symbols(circle$center_x, circle$center_y, circles = circle$radius+0.02, inches = FALSE, add = TRUE, fg = col)
#       symbols(circle$center_x, circle$center_y, circles = circle$radius-0.02, inches = FALSE, add = TRUE, fg = col)
#     }
#
#     if (valid) return(data.frame(X = circle$center_x, Y = circle$center_y, Z = circle$z, R = circle$radius, id = id))
#     else return(NULL)
#   }
#
#   clusters <- split(cl_wood@data, by = "clusterID")
#   clusters <- Filter(function(x) { nrow(x) > 20 }, clusters)
#
#   n <- length(clusters)
#   i <- 0
#   circles <- lapply(clusters, function(cl)
#   {
#     i <<- i + 1
#     if (i %% 10 == 0)  cat(sprintf("\r  Processed %d / %d", i, n))
#     fit_circle_to_seed(cl)
#   })
#
#   cat("\n")
#   circles <- Filter(Negate(is.null), circles)
#
#   circles_detected = length(circles) > 0
#
#   if (!circles_detected) {
#     warning("No circle dectected")
#   } else {
#     circles  <- do.call(rbind, circles)
#   }
#
#   if (FALSE)
#   {
#     x <- plot(cl_wood, color = 'clusterID', pal = pastel.colors(500))
#     plot(long_passages, pal = "green", add = x)
#     for (i in 1:nrow(circles))
#       add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], circles$Z[i])
#   }
#
#   logger("  Overlapping disc detection")
#
#   # Overlapping discs
#   # ------------------
#   # Pairwise distances between centers
#   circles$id = NULL
#   df = circles
#   dist_mat <- as.matrix(stats::dist(df[, c("X", "Y")], diag = TRUE, upper = TRUE))
#
#   # Define an overlap rule (e.g. centers closer than sum of radii)
#   overlap <- dist_mat < outer(df$R, df$R, "+") & dist_mat > 0
#
#   # Build adjacency graph from overlap matrix
#   n <- nrow(df)
#   visited <- logical(n)
#   group <- integer(n)
#   gid <- 0
#
#   for (i in seq_len(n)) {
#     if (!visited[i]) {
#       gid <- gid + 1
#       # breadth-first search (BFS) for connected components
#       queue <- i
#       while (length(queue) > 0) {
#         node <- queue[1]
#         queue <- queue[-1]
#         if (!visited[node]) {
#           visited[node] <- TRUE
#           group[node] <- gid
#           neighbors <- which(overlap[node, ])
#           queue <- c(queue, neighbors[!visited[neighbors]])
#         }
#       }
#     }
#   }
#
#   df$id <- group
#   df
#   circles = df
#
#   return(circles)
# }


generate_cage = function(circles, params, header)
{
  ans = generate_cage_cpp(circles,  params$path_finder$decimation)
  ans$passage = 1000
  ans$hag = 0
  cage = suppressWarnings(lidR::LAS(ans, header))
  return(cage)
}

# r_generate_cage = function(circles, params, header)
# {
#   # Convert circle to points
#   res = params$path_finder$decimation*0.75
#   circle_points_list <- lapply(1:nrow(circles), function(i)
#   {
#     generate_circle_points(circles$X[i], circles$Y[i], circles$Z[i], circles$R[i], step = res)
#   })
#   circle_points = do.call(rbind, circle_points_list)
#
#   if (FALSE) qpoints3d(circle_points)
#
#   connectors <- generate_all_connectors(circles, step_z = res)
#
#   if (FALSE) qpoints3d(connectors)
#
#   cage = rbind(circle_points, connectors)
#   cage$passage = 1000
#   cage$hag = 0
#   cage = suppressWarnings(lidR::LAS(cage, header))
#
#   if (FALSE) qpoints3d(cage@data)
#
#   return(cage)
# }
#
# # Convert to points
# generate_circle_points <- function(x, y, z, r, step = 0.1)
# {
#   # Circumference
#   circumference <- 2 * pi * r
#   # Number of points for approximately every `step` meters
#   n_points <- ceiling(circumference / step)
#   # Angles for points
#   theta <- seq(0, 2*pi, length.out = n_points + 1)[-1]  # remove last point to avoid duplicate
#   # Generate points
#   data.frame(
#     X = x + r * cos(theta),
#     Y = y + r * sin(theta),
#     Z = rep(z, n_points)
#   )
# }
#
# # --- Generate random points on a disk (XY plane) ---
# generate_disk_points <- function(x, y, z, r, n = 1000)
# {
#   theta  <- stats::runif(n, 0, 2*pi)
#   radius <- sqrt(stats::runif(n)) * r
#   data.frame(
#     X = x + radius * cos(theta),
#     Y = y + radius * sin(theta),
#     Z = rep(z, n)
#   )
# }
#
# generate_disk_radii <- function(x, y, z, r, n = 8)
# {
#   angles <- seq(0, 2*pi, length.out = n + 1)[- (n + 1)]  # remove last to avoid duplicating 0
#
#   data.frame(
#     X = x + r * cos(angles),
#     Y = y + r * sin(angles),
#     Z = rep(z, n)
#   )
# }
#
#
# # --- Generate points for one group ---
# generate_connectors <- function(df_group, step_z = 0.05)
# {
#   df_group <- df_group[order(df_group$Z), ]
#   connectors <- list()
#
#   for (i in seq_len(nrow(df_group) - 1))
#   {
#     c1 <- df_group[i, ]
#     c2 <- df_group[i + 1, ]
#
#     # vertical interpolation sequence (every step_z)
#     z_seq <- seq(c1$Z, c2$Z, by = step_z)
#     if (utils::tail(z_seq, 1) != c2$Z)
#       z_seq <- c(z_seq, c2$Z)  # include top
#
#     t_seq <- (z_seq - c1$Z) / (c2$Z - c1$Z)
#
#     # interpolate center and radius
#     x_seq <- c1$X + t_seq * (c2$X - c1$X)
#     y_seq <- c1$Y + t_seq * (c2$Y - c1$Y)
#     r_seq <- c1$R + t_seq * (c2$R - c1$R)
#
#     # add 4-radii disks
#     for (j in seq_along(z_seq))
#     {
#       connectors[[length(connectors) + 1]] <- generate_disk_radii(x_seq[j], y_seq[j], z_seq[j], r_seq[j])
#     }
#   }
#
#   connectors = do.call(rbind, connectors)
#
#   return(connectors)
# }
#
# # Wrapper for all groups
# generate_all_connectors <- function(df, step_z = 0.05)
# {
#   groups <- unique(df$id)
#   connectors <- list()
#
#   for (g in seq_along(groups))
#   {
#     df_g <- df[df$id == groups[g], ]
#     connectors[[g]] <- generate_connectors(df_g, step_z)
#   }
#
#   connectors = do.call(rbind, connectors)
#
#   if (FALSE) qpoints3d(connectors)
#
#   return(connectors)
# }


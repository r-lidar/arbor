# qsm_debug = function(tree, step = 0.2, cl_dist = 0.1, max_d = 0.1, apex = 0.0025, ...)
# {
#   #step = 0.2; cl_dist = 0.1; max_d = 0.1; apex = 0.0025
#
#   t0 <- tic() ; gc()
#
#   # Move to origin for numerical stability
#   j    <- which.min(tree$Z)
#   tx   <- tree$X[j]
#   ty   <- tree$Y[j]
#   tz   <- tree$Z[j]
#   tree <- shift(tree, tx, ty, tz)
#
#   tree <- filter_tree(tree)                        # R
#   tree <- clean_tree_butt(tree)                    # c++
#
#   qsm  <- qsm_skeleton(tree, step, cl_dist, max_d) # c++
#   qsm  <- qsm_architecture(qsm)                    # c++
#   qsm  <- qsm_smooth(qsm, niter = 1)               # c++
#   qsm  <- qsm_detect_weird_butt(qsm)               # c++
#
#   d    <- estimate_prolongation(tree, qsm)         # c++
#
#   qsm  <- qsm_prolongation(qsm, d)                 # c++
#   qsm  <- qsm_radius(qsm, tree, tip_radius = apex) # c++
#   qsm  <- qsm_volume(qsm)
#   qsm  <- shift(qsm, -tx, -ty, -tz)
#
#   data.table::setDT(qsm)
#
#   order = c("startX", "startY", "startZ", "endX", "endY", "endZ", "cyl_ID", "parent_ID", "axis_ID", "branch_order","subtree_length", "radius", "volume")
#   data.table::setcolorder(qsm, order)
#
#   qsm <- set_qsm_class(qsm)
#
#   st_crs(qsm) <- sf::st_crs(tree)
#
#   toc(t0, space = "")
#
#   return(qsm)
# }
#
# qsm_architecture <- function(qsm)
# {
#   logger("Building architecture")
#   qsm <- qsm_architecture_cpp(qsm)
#   data.table::setDT(qsm)
#   return(qsm[])
# }
#
# qsm_detect_weird_butt = function(qsm)
# {
#   axis_ID <- cyl_ID <- parent_ID <- NULL
#
#   logger("Validating butt architecture")
#
#   qsm$angle <- with(qsm,{
#     dx <- endX - startX
#     dy <- endY - startY
#     dz <- endZ - startZ
#     acos(dz / sqrt(dx^2 + dy^2 + dz^2)) * 180 / pi
#   })
#
#   main   <- qsm[axis_ID == 1]
#   angles <- main$angle
#   thresh <- 50
#   window <- 4  # number of consecutive values required below threshold
#   i      <- 1
#   while (i <= length(angles))
#   {
#     if (all(angles[i:min(i+window-1, length(angles))] < thresh)) break
#     i <- i + 1
#   }
#
#   if (i > 1)
#   {
#     logger("Detection of weird tree butt. Automatic fix triggered.", level = "WARN")
#
#     main <- qsm[axis_ID == 1]
#     rm   <- main[1:i]
#     qsm  <- qsm[!cyl_ID %in% rm$cyl_ID]
#     qsm  <- qsm_remove_disconnected_branches(qsm)
#     j <- which(qsm$axis_ID == 1)[1]
#     qsm[j, parent_ID := 0]
#   }
#
#   return(qsm)
# }
#
# clean_tree_butt = function(tree)
# {
#   logger("Cleaning tree butt")
#   data = qsm_clean_tree_butt_cpp(tree@data);
#   data.table::setDT(data)
#   tree@data = data
#   tree = lidR::las_update(tree);
#   return(tree)
#
#   # Old R code
#   # tree@data$pointID <- 1:lidR::npoints(tree)
#   # bottom <- tree[tree$Z < min(tree$Z) + 0.25]
#   # bottom <- connected_components(bottom, 0.05, 10, connectivity = 26)
#   #
#   # if (length(unique(bottom$clusterID)) > 1)
#   # {
#   #   logger("Multiple clusters at the bottom of the tree detected. Automatic cleaning triggered.", level = "WARN")
#   #
#   #   t <- table(bottom$clusterID)
#   #   i <- as.numeric(names(which.max(t)))
#   #   r <- bottom$pointID[bottom$clusterID != i]
#   #   tree <- tree[-r]
#   # }
#   #
#   # return(tree)
# }
#
# qsm_remove_disconnected_branches <- function(dt)
# {
#   axis_ID <- cyl_ID <- parent_ID <- NULL
#
#   # Ensure data.table
#   if (!data.table::is.data.table(dt)) dt <- data.table::as.data.table(dt)
#
#   # Validate required columns
#   required_cols <- c("cyl_ID", "parent_ID", "axis_ID")
#   if (!all(required_cols %in% names(dt))) {
#     stop("Input must contain columns: cyl_ID, parent_ID, axis_ID")
#   }
#
#   # Step 1: start with cylinders on the main axis
#   keep <- dt[axis_ID == 1, cyl_ID]
#   new <- keep
#
#   # Step 2: recursively find all descendants
#   repeat {
#     children <- dt[parent_ID %in% new, cyl_ID]
#     new <- setdiff(children, keep)
#     if (length(new) == 0) break
#     keep <- c(keep, new)
#   }
#
#   # Step 3: keep only connected cylinders
#   dt[cyl_ID %in% keep]
# }
#
# qsm_prolongation <- function(qsm, d, L = 0.1)
# {
#   logger("Prolongation")
#   qsm = qsm_prolongation_cpp(qsm, d, L)
#   data.table::setDT(qsm)
#   return(qsm[])
# }
#
# estimate_prolongation = function(tree, qsm)
# {
#   qsm_estimate_prolongation_cpp(tree@data, qsm)
# }
#
# qsm_radius = function(qsm, tree, tip_radius = 0.0025)
# {
#   logger("Measuring diameters")
#
#   H <- max(tree$hag)
#   R0 <- DBH_vs_H_allometry(H)/2
#
#   qsm <- qsm_conic_allometry(qsm, 2*R0, tip_radius)
#
#   if (R0 < 0.04)
#   {
#     warning("This tree is too small to be mesured. The QSM is a pure reconstruction based on allometry", call. = FALSE)
#     qsm <- qsm_conic_allometry(qsm, R0, tip_radius)
#     return(qsm)
#   }
#
#   logger("Measuring diameters")
#   qsm <- qsm_measure(tree, qsm, sarc = 180, sins = 0.2, sinl = 0.3, srmeas = 0.03)
#
#   logger("Polynomial fitting")
#   qsm <- qsm_polynomial_fitting(qsm, tip_radius)
#
#   # If we still have NAs on main axis it means interpolation failed on main
#   # axis. We have no data.
#   axis_ID <- NULL
#   main_axis <- qsm[axis_ID == 1]
#   if (anyNA(main_axis$radius))
#   {
#     warning("Not a single valid measure for this tree. The QSM is a pure reconstruction based on allometry", call. = FALSE)
#     qsm <- qsm_conic_allometry(qsm, R0, tip_radius)
#     return(qsm)
#   }
#
#   logger("Reconstruction")
#   qsm <- qsm_reconstruction(qsm, tip_radius)
#
#   return(qsm)
# }
#
# # qsm_refine_radius = function(qsm, tree)
# # {
# #   u = qsm_distances_cpp(qsm, tree@data)
# #   tree@data$r = u$radius
# #   bigcyls = filter_poi(tree, r > 0.05)
# #   bigcyls@data$cyl_ID =  as.integer(factor(bigcyls$cyl_ID)) - 1
# #   bigcyls = split(bigcyl@data, by  = "cyl_ID")
# #
# #   for (i in seq_along(bigcyls))
# #   {
# #     cyl = bigcyls[[i]][, c("X", "Y", "Z")]
# #     cyl = as.matrix(cyl)
# #     circle = fit_circle(cyl)
# #
# #     center = c(circle$center_x, circle$center_y)
# #     ffp = fit_fourier_polar(cyl, center, fill_radius = circle$radius)
# #     curve = ffp$curve
# #     curve = rbind(curve, curve[1,])
# #     pol = st_polygon(list(as.matrix(curve)))
# #     A1 = st_area(pol)
# #     A2 = pi*circle$radius^2
# #     ffp_inliers = length(ffp$inliers)
# #     ransac_inlier = length(circle$inliers)
# #     rransc = round(sqrt(A1/pi),3)
# #     rffp = round(circle$radius,3)
# #
# #     plot(cyl[, 1:2], asp = 1, main = i)
# #     symbols(circle$center_x, circle$center_y, circles = circle$radius, inches = FALSE, fg = "red", add = T, lwd = 2)
# #     symbols(circle$center_x, circle$center_y, circles = circle$radius+0.01, inches = FALSE, fg = "red", add = T, lty = 3)
# #     symbols(circle$center_x, circle$center_y, circles = circle$radius-0.01, inches = FALSE, fg = "red", add = T, lty = 3)
# #     lines(curve, col = "purple", lwd = 2)
# #     legend(
# #       "topright",
# #       bty = "n",
# #       legend = c(
# #         paste0("RANSAC radius: ", rffp),
# #         paste0("FFP equiv radius: ", rransc),
# #         paste0("RANSAC inliers: ", ransac_inlier),
# #         paste0("FFP inliers: ", ffp_inliers)
# #       )
# #     )
# #   }
# # }
#
#
# qsm_conic_allometry = function(qsm, R0, tip_radius = 0.0025)
# {
#   return(qsm_conic_allometry_cpp(qsm, R0, tip_radius))
#
#   #  old R code
#   #root = which(qsm$parent_ID == 0)
#   #w0 = qsm[["subtree_length"]][root]
#   #wi = qsm[["subtree_length"]]
#   #
#   # s = (wi / w0)
#   # s_min = min(s)
#   # s_max = max(s)
#   # r1 = tip_radius + (s - s_min) / (s_max - s_min) * (R0 - tip_radius)
#   #
#   # qsm[["radius"]] = r1
#   # return(qsm)
# }
#
# qsm_measure = function(tree, qsm, sarc = 180, sins = 0.2, sinl = 0.3, srmeas = 0.03)
# {
#   # qsm_conic_allometry must be computed first
#   qsm$theoric_radius = qsm$radius
#   theoric_radius = qsm$radius
#   qsm$radius = NULL
#
#   qsm = qsm_measure_cpp(tree@data, qsm, sarc, sins, sinl, srmeas)
#   data.table::setDT(qsm)
#   return(qsm)
# }
#
# qsm_polynomial_fitting = function(qsm, tip_radius)
# {
#   qsm = qsm_polynomial_fitting_cpp(qsm, tip_radius)
#   data.table::setDT(qsm)
#   return(qsm[])
# }
#
# qsm_reconstruction = function(qsm, tip_radius)
# {
#   qsm = qsm_reconstruction_cpp(qsm, tip_radius)
#   data.table::setDT(qsm)
#   return(qsm[])
# }
#
# # *******************************************
# # Original functions before to convert to C++
# # *******************************************
#
# qsm_measure_r = function(tree, qsm)
# {
#   # qsm_conic_allometry must be computed first
#   qsm$theoric_radius = qsm$radius
#   theoric_radius = qsm$radius
#   qsm$radius = NULL
#
#   cyl_ID <- NULL
#
#   # Assign each point to an edge of the qsm graph
#   centroid = list()
#   centroid$centroidX = (qsm$startX + qsm$endX)/2
#   centroid$centroidY = (qsm$startY + qsm$endY)/2
#   centroid$centroidZ = (qsm$startZ + qsm$endZ)/2
#   centroid = as.data.frame(centroid)
#   nearest = FNN::knnx.index(data = centroid, query = tree@data[,1:3], algorithm = "kd_tree", k = 1)
#   xyz = tree@data[, 1:3]
#   xyz$cyl_ID = qsm$cyl_ID[nearest]
#
#   if (FALSE)
#   {
#     x = plot_qsm(qsm, cylinder = FALSE, add = c(0,0))
#     cc = centroid
#     rgl::points3d(cc, size = 3)
#     rgl::points3d(xyz, col = lidR:::set.colors(xyz$cyl_ID, pal = lidR::random.colors(nrow(centroid))))
#     lidR:::.pan3d(2)
#   }
#
#   # Split data by cylinder (graph edges)
#   data.table::setkey(qsm, cyl_ID)
#   data.table::setkey(xyz, cyl_ID)
#
#   r2 = rep(NA_real_, nrow(qsm))
#   r_rescue = rep(NA_real_, nrow(qsm))
#
#   subs = split(xyz, by = "cyl_ID",  keep.by = FALSE)
#   starts = qsm[, c(1:3, 7)]
#   ends = qsm[, c(4:6, 7)]
#   starts = split(starts, by = "cyl_ID",  keep.by = FALSE)
#   ends = split(ends, by = "cyl_ID",  keep.by = FALSE)
#   theoric_radii = qsm[, c("cyl_ID", "theoric_radius")]
#   theoric_radii = split(theoric_radii, by = "cyl_ID",  keep.by = FALSE)
#
#   # Actual measurement for each cylindre using RANSAC if it works
#   is.valid.circle = function(radius, angle_range, pinliner, pinside)
#   {
#     if (radius < 0.03)  return(FALSE)
#     if (pinside > 20) return(FALSE)
#     return(angle_range > 120 & pinliner > 30)
#   }
#
#   for (i in 1:nrow(qsm))
#   {
#     id = qsm$cyl_ID[i]
#
#     if (id < 1) next # skip prolongations
#
#     id = as.character(id)
#
#     sub = subs[[id]]
#
#     if (is.null(sub)) next
#
#     theoric_radius = as.numeric(theoric_radii[[id]])
#
#     # No need to try to measure anything if the theory is < 4 cm
#     if (theoric_radius < 0.04) next
#
#     sub = as.matrix(sub)
#
#     start = as.numeric(starts[[id]])
#     end = as.numeric(ends[[id]])
#
#     if (FALSE)
#     {
#       rgl::points3d(sub)
#       rgl::axes3d()
#       rgl::segments3d(rbind(start, end), lwd = 3)
#     }
#
#     M = compute_rotation_matrix(start, end)
#     sub = sub %*% t(M)
#
#     if (nrow(sub) < 50) next
#
#     r = ransac_circle(sub, num_iterations = 100, inlier_threshold = 0.02)
#
#     valid = is.valid.circle(r$radius, r$covered_arc_degree, r$percentage_inlier*100, r$percentage_inside*100)
#
#     if (valid)
#     {
#       r2[i] = r$radius
#
#       if (FALSE)
#       {
#         col = if(valid) "darkgreen" else "red"
#         title = paste0("ID = ",  i, " | arc = ", round(r$covered_arc_degree, 1),  " | inliner = ", round(r$percentage_inlier*100), "% insider = ", round(r$percentage_inside*100), "% | R = ", round(r$radius, 2))
#         plot(sub, asp = 1, pch = 19, cex = 0.5, main = title)
#         points(r$center_x, r$center_y, pch = 4)
#         symbols(r$center_x, r$center_y, circles = r$radius, add = TRUE, fg = col, inches = FALSE, lwd = 2)
#         symbols(r$center_x, r$center_y, circles = r$radius+0.01, add = TRUE, fg = col, inches = FALSE, lty = 3, lwd = 2)
#         symbols(r$center_x, r$center_y, circles = r$radius-0.01, add = TRUE, fg = col, inches = FALSE, lty = 3, lwd = 2)
#       }
#     }
#   }
#
#   qsm[["radius"]] = data.table::copy(r2)
#
#   if (FALSE)
#   {
#     x = plot_qsm(qsm, cylinder = TRUE)
#     plot(tree, add = x, size = 2, pal = "brown")
#   }
#
#   return(qsm)
# }
#
# qsm_polynomial_fitting_r = function(qsm, tip_radius)
# {
#   radius <- axis_ID <- NULL
#
#   # We now have some cynlindred measure with ransac
#   # Polynomial interpolation for each axis with enought measurements
#   for (i in sort(unique(qsm$axis_ID)))
#   {
#     axe = qsm[axis_ID == i]
#     if (sum(!is.na(axe$radius)) <= 6) next
#
#     mes = axe[!is.na(radius)]
#     mod <- stats::nls(
#       radius ~ tip_radius + a * subtree_length + b * subtree_length^2,
#       data = mes,
#       start = list(a = 0.01, b = 0.001),
#       algorithm = "port",
#       lower = c(-Inf, -Inf),
#       upper = c(Inf, Inf)   # b > 0
#     )
#
#     r = stats::predict(mod, axe)
#     qsm[axis_ID == i, radius := r]
#
#     if (FALSE)
#     {
#       # plot(axe$theoric_radius, max(axe$subtree_length) -axe$subtree_length,
#       #      asp = 0.01,
#       #      main = NULL,
#       #      pch  = 19,
#       #      cex = 0.5,
#       #      col = "gray",
#       #      xlab = "Radius (m)",
#       #      ylab = "Distance to root (m)",
#       #      xlim = c(0, 0.2))
#       plot(axe$radius, max(axe$subtree_length) -axe$subtree_length,
#            pch = 19, cex = 0.6,
#            xlab = "Radius (m)",
#            ylab = "Distance to root (m)",
#            xlim = c(0, max(r, na.omit(axe$radius))))
#       points(r, max(axe$subtree_length) -axe$subtree_length, col = "blue", pch = 19, cex = 0.5)
#
#       legend("topright",
#              legend = c("Theoric radius", "Measured radius", "Interpolation"),
#              col = c("gray", "black", "blue"),
#              pch = 19,
#              pt.cex = c(0.5, 0.6, 0.5))
#     }
#   }
#
#   return(qsm[])
# }
#
# qsm_reconstruction_r = function(qsm, tip_radius)
# {
#   axis_ID <- radius <- branch_order <- NULL
#
#   # We have some cylinder measured and interpolated. A large portion of the tree is missing.
#   # We loop on each axes by branch order. If we have some NAs we compare to the theory.
#   # We ensure the theory is not bigger than the parent branch. If the theory is bigger
#   # We rescale to be the size of the parent
#   # And if we have enough measurement we scale the theory to something realistic
#
#   # Loop by branch order
#   for (i in sort(unique(qsm$branch_order)))
#   {
#     if (i == 1) next
#
#     order = qsm[branch_order == i]
#     axis_IDs = unique(order$axis_ID)
#
#     for (j in axis_IDs)
#     {
#       axe = order[axis_ID == j]
#
#       # Find parent radius and subtree length
#       parent = axe$parent_ID[1]
#       parent_id = which(qsm$cyl_ID == parent)
#       if (length(parent_id) == 0) stop("Internal error. No parent.")
#
#       # Compute the theoretical radius
#       r0 = qsm$radius[parent_id]
#       w0 = qsm$subtree_length[parent_id]
#       wi = axe$subtree_length
#       s  = (wi / w0)^1.1
#       r  = tip_radius + s * (r0 - tip_radius)
#
#       # We have only NAs. We can only use the theory
#       if (all(is.na(axe$radius)))
#       {
#
#       }
#       # We have some measures. Use a comparison between measure and theory to rescale theory
#       else if(anyNA(axe$radius))
#       {
#         r_mes = axe$radius
#
#         # We need at least 3 measurements otherwise it is easy to have
#         # false ransac fitting. Other wise use radius theoretic
#         if (sum(!is.na(r_mes)) >= 3)
#         {
#           # Look at the average ratio between the theory in and the measures.
#           # If measures tends to be bigger increase theory. Or opposite
#           ratio = stats::median(r_mes/r, na.rm = TRUE)
#           r = tip_radius + s * (r0*ratio - tip_radius)
#         }
#       }
#
#       qsm[branch_order == i & axis_ID == j, radius := r]
#     }
#   }
#
#   if (FALSE)
#   {
#     x = plot_qsm(qsm, cylinder = TRUE)
#     plot(tree, add =x, size = 2, pal = "brown")
#   }
#
#   return(qsm[])
# }
#
# qsm_skeleton = function(tree, step = .2, cl_dist = 0.1, max_d = 0.3)
# {
#   data <- qsm_layers(tree, step)
#   data <- qsm_clusters(data, cl_dist)
#   skel <- qsm_nodes(data, max_d)
#   skel <- qsm_topology(skel)
#   return(skel)
# }
#
# qsm_layers = function(tree, step)
# {
#   logger("Computing layers")
#   data = qsm_layers_cpp(tree@data, step)
#
#   res = list()
#   res[["X"]] = tree@data[["X"]]
#   res[["Y"]] = tree@data[["Y"]]
#   res[["Z"]] = tree@data[["Z"]]
#   res[["iter"]] = data[["iter"]]
#   res[["dist"]] = data[["dist"]]
#   data.table::setDT(res)
#   return(res[])
# }
#
# qsm_clusters = function(data, cl_dist)
# {
#   res = qsm_cluster_cpp(data, cl_dist)
#   data[["cluster"]] = res[["cluster"]]
#   data[["radius"]] = res[["radius"]]
#   return(data)
# }
#
# # @importFrom data.table :=
# # qsm_clusters = function(data, cl_dist)
# # {
# #   iter <- cluster <- radius <- X <- Y <- Z <- NULL
# #
# #   logger("Clustering layers")
# #
# #   first = TRUE
# #   for (i in sort(unique(data$iter)))
# #   {
# #     if (first)
# #     {
# #       # at the first iteration the clustering distance (cl_dist) is the radius of the first layer
# #       # only one cluster at iteration 1 (the tree base)
# #       data[iter == i, cluster := 1]
# #       # compute the radius
# #       data[iter == i, radius := sqrt((X-mean(X))^2 + (Y-mean(Y))^2 + (Z-mean(Z))^2  )]
# #       # the average radius is the first clustering distance
# #       cl_d = mean(data$radius[data$iter == i])
# #       first = FALSE
# #
# #       if (FALSE)
# #       {
# #         tmp = data[iter == i]
# #         plot(tmp[, .(X,Y)], asp = 1, cex = 0.25)
# #         points(mean(tmp$X), mean(tmp$Y), col = "red")
# #       }
# #     }
# #     else
# #     {
# #       # at subsequent iterations the clustering distance is the average radius
# #       # of clustered objects in the curent layer
# #
# #       # keep points in the curent layer
# #       in_iter = which(data$iter==i)
# #
# #       # if there are more than one point in the layer -> cluster it
# #       if (length(in_iter) >= 2)
# #       {
# #         # clustering
# #         cl = dbscan::dbscan(data[in_iter,1:3], eps = cl_dist, minPts = 1)
# #         #cl = lidR::connected_components(temp, res = cl_dist, min_pts = 1, connectivity = 26)
# #         cl_id = cl$cluster
# #         data[in_iter, cluster := cl_id]
# #         data[in_iter, radius := sqrt((X-mean(X))^2 + (Y-mean(Y))^2 + (Z-mean(Z))^2  ), by = cluster]  # compute the radius for each cluster
# #         cl_d = mean(data$radius[data$iter == i])/2           # the new clustering distance is the average radius of all clusters
# #
# #         if (FALSE)
# #         {
# #           tmp = data[in_iter]
# #           plot(tmp[, 1:2], asp = 1, col = cl_id, cex = 0.25)
# #         }
# #
# #         # if the computed clustering distance is too small -> replace by the user defined minimum distance
# #         #if(cl_d < cl_dist) cl_d = cl_dist
# #       }
# #       else
# #       {
# #         # if there is only one point, there is only one cluster
# #         data[in_iter,cluster := 1]
# #       }
# #     }
# #   }
# #
# #   if (FALSE)
# #   {
# #     col = lidR:::set.colors(data$cluster+data$iter, pastel.colors(50000))
# #     rgl::points3d(data$X, data$Y, data$Z, col = col)
# #     col = lidR:::set.colors(data$radius, height.colors(500))
# #     rgl::points3d(data$X, data$Y, data$Z, col = col)
# #   }
# #
# #   return(data[])
# # }
#
# qsm_nodes = function(data, max_d)
# {
#   logger("Building nodes")
#   skel = cpp_build_skeleton(data, max_d)
#   data.table::setDT(skel)
#
#   if (FALSE)
#   {
#     x = plot(tree, bg = "white")
#     plot_qsm(skel, add = x)
#   }
#
#   return(skel)
# }
#
# qsm_radius = function(qsm, tree, tip_radius = 0.0025)
# {
#   logger("Measuring diameters")
#
#   H <- max(tree$hag)
#   R0 <- DBH_vs_H_allometry(H)/2
#
#   qsm <- qsm_conic_allometry(qsm, 2*R0, tip_radius)
#
#   if (R0 < 0.04)
#   {
#     warning("This tree is too small to be mesured. The QSM is a pure reconstruction based on allometry", call. = FALSE)
#     qsm <- qsm_conic_allometry(qsm, R0, tip_radius)
#     return(qsm)
#   }
#
#   logger("Measuring diameters")
#   qsm <- qsm_measure(tree, qsm, sarc = 180, sins = 0.2, sinl = 0.3, srmeas = 0.03)
#
#   logger("Polynomial fitting")
#   qsm <- qsm_polynomial_fitting(qsm, tip_radius)
#
#   # If we still have NAs on main axis it means interpolation failed on main
#   # axis. We have no data.
#   axis_ID <- NULL
#   main_axis <- qsm[axis_ID == 1]
#   if (anyNA(main_axis$radius))
#   {
#     warning("Not a single valid measure for this tree. The QSM is a pure reconstruction based on allometry", call. = FALSE)
#     qsm <- qsm_conic_allometry(qsm, R0, tip_radius)
#     return(qsm)
#   }
#
#   logger("Reconstruction")
#   qsm <- qsm_reconstruction(qsm, tip_radius)
#
#   return(qsm)
# }
#
# # qsm_refine_radius = function(qsm, tree)
# # {
# #   u = qsm_distances_cpp(qsm, tree@data)
# #   tree@data$r = u$radius
# #   bigcyls = filter_poi(tree, r > 0.05)
# #   bigcyls@data$cyl_ID =  as.integer(factor(bigcyls$cyl_ID)) - 1
# #   bigcyls = split(bigcyl@data, by  = "cyl_ID")
# #
# #   for (i in seq_along(bigcyls))
# #   {
# #     cyl = bigcyls[[i]][, c("X", "Y", "Z")]
# #     cyl = as.matrix(cyl)
# #     circle = fit_circle(cyl)
# #
# #     center = c(circle$center_x, circle$center_y)
# #     ffp = fit_fourier_polar(cyl, center, fill_radius = circle$radius)
# #     curve = ffp$curve
# #     curve = rbind(curve, curve[1,])
# #     pol = st_polygon(list(as.matrix(curve)))
# #     A1 = st_area(pol)
# #     A2 = pi*circle$radius^2
# #     ffp_inliers = length(ffp$inliers)
# #     ransac_inlier = length(circle$inliers)
# #     rransc = round(sqrt(A1/pi),3)
# #     rffp = round(circle$radius,3)
# #
# #     plot(cyl[, 1:2], asp = 1, main = i)
# #     symbols(circle$center_x, circle$center_y, circles = circle$radius, inches = FALSE, fg = "red", add = T, lwd = 2)
# #     symbols(circle$center_x, circle$center_y, circles = circle$radius+0.01, inches = FALSE, fg = "red", add = T, lty = 3)
# #     symbols(circle$center_x, circle$center_y, circles = circle$radius-0.01, inches = FALSE, fg = "red", add = T, lty = 3)
# #     lines(curve, col = "purple", lwd = 2)
# #     legend(
# #       "topright",
# #       bty = "n",
# #       legend = c(
# #         paste0("RANSAC radius: ", rffp),
# #         paste0("FFP equiv radius: ", rransc),
# #         paste0("RANSAC inliers: ", ransac_inlier),
# #         paste0("FFP inliers: ", ffp_inliers)
# #       )
# #     )
# #   }
# # }
#
#
# qsm_conic_allometry = function(qsm, R0, tip_radius = 0.0025)
# {
#   return(qsm_conic_allometry_cpp(qsm, R0, tip_radius))
#
#   #  old R code
#   #root = which(qsm$parent_ID == 0)
#   #w0 = qsm[["subtree_length"]][root]
#   #wi = qsm[["subtree_length"]]
#   #
#   # s = (wi / w0)
#   # s_min = min(s)
#   # s_max = max(s)
#   # r1 = tip_radius + (s - s_min) / (s_max - s_min) * (R0 - tip_radius)
#   #
#   # qsm[["radius"]] = r1
#   # return(qsm)
# }
#
# qsm_measure = function(tree, qsm, sarc = 180, sins = 0.2, sinl = 0.3, srmeas = 0.03)
# {
#   # qsm_conic_allometry must be computed first
#   qsm$theoric_radius = qsm$radius
#   theoric_radius = qsm$radius
#   qsm$radius = NULL
#
#   qsm = qsm_measure_cpp(tree@data, qsm, sarc, sins, sinl, srmeas)
#   data.table::setDT(qsm)
#   return(qsm)
# }
#
# qsm_polynomial_fitting = function(qsm, tip_radius)
# {
#   qsm = qsm_polynomial_fitting_cpp(qsm, tip_radius)
#   data.table::setDT(qsm)
#   return(qsm[])
# }
#
# qsm_reconstruction = function(qsm, tip_radius)
# {
#   qsm = qsm_reconstruction_cpp(qsm, tip_radius)
#   data.table::setDT(qsm)
#   return(qsm[])
# }
#
# # *******************************************
# # Original functions before to convert to C++
# # *******************************************
#
# qsm_measure_r = function(tree, qsm)
# {
#   # qsm_conic_allometry must be computed first
#   qsm$theoric_radius = qsm$radius
#   theoric_radius = qsm$radius
#   qsm$radius = NULL
#
#   cyl_ID <- NULL
#
#   # Assign each point to an edge of the qsm graph
#   centroid = list()
#   centroid$centroidX = (qsm$startX + qsm$endX)/2
#   centroid$centroidY = (qsm$startY + qsm$endY)/2
#   centroid$centroidZ = (qsm$startZ + qsm$endZ)/2
#   centroid = as.data.frame(centroid)
#   nearest = FNN::knnx.index(data = centroid, query = tree@data[,1:3], algorithm = "kd_tree", k = 1)
#   xyz = tree@data[, 1:3]
#   xyz$cyl_ID = qsm$cyl_ID[nearest]
#
#   if (FALSE)
#   {
#     x = plot_qsm(qsm, cylinder = FALSE, add = c(0,0))
#     cc = centroid
#     rgl::points3d(cc, size = 3)
#     rgl::points3d(xyz, col = lidR:::set.colors(xyz$cyl_ID, pal = lidR::random.colors(nrow(centroid))))
#     lidR:::.pan3d(2)
#   }
#
#   # Split data by cylinder (graph edges)
#   data.table::setkey(qsm, cyl_ID)
#   data.table::setkey(xyz, cyl_ID)
#
#   r2 = rep(NA_real_, nrow(qsm))
#   r_rescue = rep(NA_real_, nrow(qsm))
#
#   subs = split(xyz, by = "cyl_ID",  keep.by = FALSE)
#   starts = qsm[, c(1:3, 7)]
#   ends = qsm[, c(4:6, 7)]
#   starts = split(starts, by = "cyl_ID",  keep.by = FALSE)
#   ends = split(ends, by = "cyl_ID",  keep.by = FALSE)
#   theoric_radii = qsm[, c("cyl_ID", "theoric_radius")]
#   theoric_radii = split(theoric_radii, by = "cyl_ID",  keep.by = FALSE)
#
#   # Actual measurement for each cylindre using RANSAC if it works
#   is.valid.circle = function(radius, angle_range, pinliner, pinside)
#   {
#     if (radius < 0.03)  return(FALSE)
#     if (pinside > 20) return(FALSE)
#     return(angle_range > 120 & pinliner > 30)
#   }
#
#   for (i in 1:nrow(qsm))
#   {
#     id = qsm$cyl_ID[i]
#
#     if (id < 1) next # skip prolongations
#
#     id = as.character(id)
#
#     sub = subs[[id]]
#
#     if (is.null(sub)) next
#
#     theoric_radius = as.numeric(theoric_radii[[id]])
#
#     # No need to try to measure anything if the theory is < 4 cm
#     if (theoric_radius < 0.04) next
#
#     sub = as.matrix(sub)
#
#     start = as.numeric(starts[[id]])
#     end = as.numeric(ends[[id]])
#
#     if (FALSE)
#     {
#       rgl::points3d(sub)
#       rgl::axes3d()
#       rgl::segments3d(rbind(start, end), lwd = 3)
#     }
#
#     M = compute_rotation_matrix(start, end)
#     sub = sub %*% t(M)
#
#     if (nrow(sub) < 50) next
#
#     r = ransac_circle(sub, num_iterations = 100, inlier_threshold = 0.02)
#
#     valid = is.valid.circle(r$radius, r$covered_arc_degree, r$percentage_inlier*100, r$percentage_inside*100)
#
#     if (valid)
#     {
#       r2[i] = r$radius
#
#       if (FALSE)
#       {
#         col = if(valid) "darkgreen" else "red"
#         title = paste0("ID = ",  i, " | arc = ", round(r$covered_arc_degree, 1),  " | inliner = ", round(r$percentage_inlier*100), "% insider = ", round(r$percentage_inside*100), "% | R = ", round(r$radius, 2))
#         plot(sub, asp = 1, pch = 19, cex = 0.5, main = title)
#         points(r$center_x, r$center_y, pch = 4)
#         symbols(r$center_x, r$center_y, circles = r$radius, add = TRUE, fg = col, inches = FALSE, lwd = 2)
#         symbols(r$center_x, r$center_y, circles = r$radius+0.01, add = TRUE, fg = col, inches = FALSE, lty = 3, lwd = 2)
#         symbols(r$center_x, r$center_y, circles = r$radius-0.01, add = TRUE, fg = col, inches = FALSE, lty = 3, lwd = 2)
#       }
#     }
#   }
#
#   qsm[["radius"]] = data.table::copy(r2)
#
#   if (FALSE)
#   {
#     x = plot_qsm(qsm, cylinder = TRUE)
#     plot(tree, add = x, size = 2, pal = "brown")
#   }
#
#   return(qsm)
# }
#
# qsm_polynomial_fitting_r = function(qsm, tip_radius)
# {
#   radius <- axis_ID <- NULL
#
#   # We now have some cynlindred measure with ransac
#   # Polynomial interpolation for each axis with enought measurements
#   for (i in sort(unique(qsm$axis_ID)))
#   {
#     axe = qsm[axis_ID == i]
#     if (sum(!is.na(axe$radius)) <= 6) next
#
#     mes = axe[!is.na(radius)]
#     mod <- stats::nls(
#       radius ~ tip_radius + a * subtree_length + b * subtree_length^2,
#       data = mes,
#       start = list(a = 0.01, b = 0.001),
#       algorithm = "port",
#       lower = c(-Inf, -Inf),
#       upper = c(Inf, Inf)   # b > 0
#     )
#
#     r = stats::predict(mod, axe)
#     qsm[axis_ID == i, radius := r]
#
#     if (FALSE)
#     {
#       # plot(axe$theoric_radius, max(axe$subtree_length) -axe$subtree_length,
#       #      asp = 0.01,
#       #      main = NULL,
#       #      pch  = 19,
#       #      cex = 0.5,
#       #      col = "gray",
#       #      xlab = "Radius (m)",
#       #      ylab = "Distance to root (m)",
#       #      xlim = c(0, 0.2))
#       plot(axe$radius, max(axe$subtree_length) -axe$subtree_length,
#            pch = 19, cex = 0.6,
#            xlab = "Radius (m)",
#            ylab = "Distance to root (m)",
#            xlim = c(0, max(r, na.omit(axe$radius))))
#       points(r, max(axe$subtree_length) -axe$subtree_length, col = "blue", pch = 19, cex = 0.5)
#
#       legend("topright",
#              legend = c("Theoric radius", "Measured radius", "Interpolation"),
#              col = c("gray", "black", "blue"),
#              pch = 19,
#              pt.cex = c(0.5, 0.6, 0.5))
#     }
#   }
#
#   return(qsm[])
# }
#
# qsm_reconstruction_r = function(qsm, tip_radius)
# {
#   axis_ID <- radius <- branch_order <- NULL
#
#   # We have some cylinder measured and interpolated. A large portion of the tree is missing.
#   # We loop on each axes by branch order. If we have some NAs we compare to the theory.
#   # We ensure the theory is not bigger than the parent branch. If the theory is bigger
#   # We rescale to be the size of the parent
#   # And if we have enough measurement we scale the theory to something realistic
#
#   # Loop by branch order
#   for (i in sort(unique(qsm$branch_order)))
#   {
#     if (i == 1) next
#
#     order = qsm[branch_order == i]
#     axis_IDs = unique(order$axis_ID)
#
#     for (j in axis_IDs)
#     {
#       axe = order[axis_ID == j]
#
#       # Find parent radius and subtree length
#       parent = axe$parent_ID[1]
#       parent_id = which(qsm$cyl_ID == parent)
#       if (length(parent_id) == 0) stop("Internal error. No parent.")
#
#       # Compute the theoretical radius
#       r0 = qsm$radius[parent_id]
#       w0 = qsm$subtree_length[parent_id]
#       wi = axe$subtree_length
#       s  = (wi / w0)^1.1
#       r  = tip_radius + s * (r0 - tip_radius)
#
#       # We have only NAs. We can only use the theory
#       if (all(is.na(axe$radius)))
#       {
#
#       }
#       # We have some measures. Use a comparison between measure and theory to rescale theory
#       else if(anyNA(axe$radius))
#       {
#         r_mes = axe$radius
#
#         # We need at least 3 measurements otherwise it is easy to have
#         # false ransac fitting. Other wise use radius theoretic
#         if (sum(!is.na(r_mes)) >= 3)
#         {
#           # Look at the average ratio between the theory in and the measures.
#           # If measures tends to be bigger increase theory. Or opposite
#           ratio = stats::median(r_mes/r, na.rm = TRUE)
#           r = tip_radius + s * (r0*ratio - tip_radius)
#         }
#       }
#
#       qsm[branch_order == i & axis_ID == j, radius := r]
#     }
#   }
#
#   if (FALSE)
#   {
#     x = plot_qsm(qsm, cylinder = TRUE)
#     plot(tree, add =x, size = 2, pal = "brown")
#   }
#
#   return(qsm[])
# }
#
# qsm_smooth = function(qsm, niter = 2)
# {
#   ans = qsm_smooth_cpp(qsm, niter = niter)
#   data.table::setDT(ans)
#   return(ans)
# }
#
# qsm_topology = function(skeleton)
# {
#   logger("Computing qsm topology")
#
#   skeleton$cyl_ID = 1:nrow(skeleton)
#   skeleton$parent_ID = 0
#   skeleton = qsm_topology_cpp(skeleton)
#
#   n_root = sum(skeleton$parent_ID == 0)
#
#   # No root is a bug
#   if (n_root == 0)
#     stop("Internal error: no root found")
#
#   # Two roots are rare but possible if forking at root
#   # Add a 1 mm cylinder to force one root
#   if (n_root > 1)
#   {
#     xyz = skeleton[which(skeleton$parent_ID == 0)[1],]
#     xyz$endX = xyz$startX
#     xyz$endY = xyz$startY
#     xyz$endZ = xyz$startZ
#     xyz$startZ = xyz$startZ - 0.001
#     skeleton = rbind(xyz, skeleton)
#     skeleton$cyl_ID = 1:nrow(skeleton)
#     skeleton$parent_ID = 0
#     skeleton = qsm_topology_cpp(skeleton)
#   }
#
#   data.table::setDT(skeleton)
#
#   return(skeleton)
# }
#

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

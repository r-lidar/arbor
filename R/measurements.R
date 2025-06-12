#' Measure the height of each tree
#'
#' @param trees LAS object with individual tree segmented.
#' @return data.frame
#' @export
measure_heights = function(trees)
{
  attributes = names(trees)
  stopifnot("treeID" %in% attributes)

  . <- .I <- treeID <- X <- Y <- Z <- NULL

  ans = trees@data[, list(Height = max(Z), id = .I[which.max(Z)]), by = treeID]
  tops = trees[ans$id]
  tops =  tops@data[, .(X,Y,Z)]
  tops$treeID = ans$treeID
  tops
}

#' Measure the diameter of each tree
#'
#' Fits a circle to each tree. This is not a DBH measurement. DBH is ill-defined in too many cases for
#' example when the tree forks lower than 1.3 m, when the trees are too small to robustly measure a DBH
#' and in many other edges cases actually. The methods actually fits numerous circle on numerous slices
#' of the trees and retains those that best fit. Ultimately it returns the circles that have the
#' largest diameters for each tree and these circles are likely to be at the bottom of the trees but
#' they may be somewhere else if it fits better.
#'
#' @param ... unused. Only useful to separate important parameters from those that should only be changed
#' in very specific contexts.
#' @param max_height Make slices for the bottom of the tree (i.e., 50 cm if you slice at 50 cm) up to
#' +max_height above that point (i.e., 1.5 m in this example).
#' @param min_slice_thickness,max_slice_thickness The algorithm fits circles on many slices with various
#' thicknesses. These control the thickness limits.
#' @param debug Boolean. Some debugging display.
#' @param trees LAS object with individual tree segmented.
#' @export
measure_diameters = function(trees, ..., max_height = 2, min_slice_thickness = 0.1, max_slice_thickness = 0.6, debug = FALSE)
{
  .<- X <- Y<- Z <- treeID <- foliage <- hag <- z <- radius <- bottom <- NULL

  attributes = names(trees)
  stopifnot("foliage" %in% attributes)
  stopifnot("hag" %in% attributes)
  stopifnot("treeID" %in% attributes)

  # Generate multiple slices
  ranges = expand.grid(bottom = seq(0.1, max_height, 0.1), top = seq(0.1, max_height, 0.1))
  ranges = ranges[ranges$top - ranges$bottom >= min_slice_thickness,]
  ranges = ranges[ranges$top - ranges$bottom <= max_slice_thickness,]
  data.table::setDT(ranges)
  data.table::setorder(ranges, bottom)

  cat("RANSAC fitting on", nrow(ranges), "slices per trees\n")

  extensions = list()
  for (id in unique(trees$treeID))
  {
    cat("Tree", id, ": ")
    tt = lidR::filter_poi(trees, treeID == id, foliage == FALSE)
    min_hag = min(tt$hag)
    tt = lidR::filter_poi(tt, hag < min_hag + max_height)
    tt = decimate_points(tt, random_per_voxel(0.02))


    #xyz = sf::st_coordinates(tt)
    #pca <- prcomp(xyz, center = TRUE, scale. = FALSE)
    #main_axis <- pca$rotation[, 1]  # First principal component

    # Compute the rotation matrix
    #rotation_matrix <- align_to_z(main_axis)

    # Apply the rotation to the point cloud
    #rotated_xyz = xyz %*% t(rotation_matrix)
    #rotated_xyz <- scale(rotated_xyz, center = TRUE, scale = FALSE)
    #ttt <- as.data.frame(rotated_xyz)
    #names(ttt) = c("X", "Y", "Z")
    #ttt$clusterID = tt$clusterID
    #ttt = LAS(ttt)


    # Plot the point cloud
    #rgl::plot3d(centered_xyz, col = "blue", size = 2)
    #rgl::points3d(rotated_xyz, col = "red", size = 2)
    ##rgl::arrow3d(p0 = c(0,0,0), p1 =  main_axis, type = "lines",  col = "red", length = 2)
    sum_valid = 0
    circles = apply(ranges, 1, function(x)
    {
      if (sum_valid > 10) return(NULL)

      bottom = tt@data[hag  >= x[1] & hag <= x[2], .(X,Y,Z)]
      bottom = as.matrix(bottom)
      if (nrow(bottom) < 10) return(NULL)

      circle = ransac_circle_cpp(bottom, early_exit = 0.7)
      valid = is.valid.circle(circle$radius, circle$covered_arc_degree, circle$percentage_inlier*100, circle$percentage_inside*100)

      if (valid) sum_valid <<- sum_valid+1

      circle$valid = valid
      circle$z = circle$z + diff(x)/2
      circle$inliers = NULL

      if (debug)
      {
        col = ifelse(valid, "darkgreen", "red")
        plot(bottom, asp = 1, main = id, pch = 19, cex = 0.25)
        graphics::symbols(circle$center_x, circle$center_y, circles = circle$radius,  add = TRUE, fg = col, inches = FALSE)
        graphics::symbols(circle$center_x, circle$center_y, circles = circle$radius+0.01,  add = TRUE, fg = col, lty=3, inches = FALSE)
        graphics::symbols(circle$center_x, circle$center_y, circles = circle$radius-0.01,  add = TRUE, fg = col, lty= 3, inches = FALSE)
        graphics::mtext(paste0("Radius = ", round(circle$radius*100, 1), " cm | inliner = ", round(circle$percentage_inlier*100), "% | sector ", circle$covered_arc_degree, " deg | HAG = [", x[1], ",", x[2], "]"))
      }

      as.data.frame(circle)
    })

    if (is.null(circles))
    {
      cat("no valid circle detected\n")
      next
    }

    circles = data.table::rbindlist(circles)

    n = sum(circles$valid)
    if (n == 0)
    {
      cat("no valid circle detected\n")
      next
    }

    circles = circles[circles$valid,]
    data.table::setorder(circles, -radius)
    n = min(c(5, n))

    if (n == 0) next

    circle = list(center_x = circles$center_x[1],
                  center_y = circles$center_y[1],
                  center_z = circles$z[1],
                  radius = circles$radius[1])

    extensions[[as.character(id)]] = circle
    cat("radius =", round(circle$radius*100, 1), "cm")

    cat("\n")
    #x = plot(tt)
    #plot(LAS(extension), add = x)
  }

  treeID = as.numeric(names(extensions))
  extensions = data.table::rbindlist(extensions)
  extensions$treeID = treeID

  return(extensions)
}




measure_heights = function(trees)
{
  ans = trees@data[, list(Height = max(Z), id = .I[which.max(Z)]), by = treeID]
  tops = trees[ans$id]
  tops =  tops@data[, .(X,Y,Z)]
  tops$treeID = ans$treeID
  tops
}

measure_diameters = function(trees, ..., max_height = 2, min_slice_thickness = 0.1, max_slice_thickness = 0.6, debug = FALSE)
{
  treeID <- foliage <- hag <- NULL

  # Generate multiple slices
  ranges = expand.grid(bottom = seq(0.1, max_height, 0.05), top = seq(0.1, max_height, 0.05))
  ranges = ranges[ranges$top - ranges$bottom >= min_slice_thickness,]
  ranges = ranges[ranges$top - ranges$bottom <= max_slice_thickness,]

  cat("RANSAC fitting on", nrow(ranges), "slices per trees\n")

  extensions = list()
  for (id in unique(trees$treeID))
  {
    cat("Tree", id, ": ")
    tt = lidR::filter_poi(trees, treeID == id, foliage == FALSE)
    min_hag = min(tt$hag)
    tt = lidR::filter_poi(tt, hag < min_hag + max_height)

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

    circles = apply(ranges, 1, function(x)
    {
      bottom = tt@data[hag  >= x[1] & hag <= x[2], .(X,Y, Z)]
      bottom = as.matrix(bottom)
      if (nrow(bottom) < 10) return(NULL)

      circle = .fit_circle(bottom)
      valid = is.valid.circle(circle$radius, circle$covered_arc_degree, circle$percentage_inlier*100)

      circle$valid = valid
      circle$z = circle$z + diff(x)/2


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
    data.table::setorder(circles, z, -radius)
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

.fit_circle_on_3_points <- function(points_subset)
{
  stopifnot(nrow(points_subset) == 3L, ncol(points_subset) > 2L)

  # Extract the coordinates
  x1 <- points_subset[1, 1]
  y1 <- points_subset[1, 2]
  x2 <- points_subset[2, 1]
  y2 <- points_subset[2, 2]
  x3 <- points_subset[3, 1]
  y3 <- points_subset[3, 2]

  # Calculate the coefficients for the linear system
  A <- 2 * (x2 - x1)
  B <- 2 * (y2 - y1)
  C <- x2^2 + y2^2 - x1^2 - y1^2
  D <- 2 * (x3 - x1)
  E <- 2 * (y3 - y1)
  G <- x3^2 + y3^2 - x1^2 - y1^2

  # Solve for a and b using Cramer's rule
  denominator <- A * E - B * D
  if (denominator == 0)
  {
    return(c(0, 0, 0))
  }
  a <- (C * E - B * G) / denominator
  b <- (A * G - C * D) / denominator

  # Calculate the radius
  r <- sqrt((x1 - a)^2 + (y1 - b)^2)

  # Return the center and radius
  return(c(a, b, r))
}

.fit_circle <- function(points, num_iterations = 100, inlier_threshold = 0.01)
{
  best_circle <- NULL
  max_inliers <- 0

  if (is(points, "LAS")) points = sf::st_coordinates(points)

  stopifnot(is.matrix(points), ncol(points) == 3L, nrow(points) > 3L)

  z = points[, 3]

  for (i in 1:num_iterations)
  {
    # Randomly sample points
    sample_indices <- sample(1:nrow(points), 3L)
    points_subset <- points[sample_indices, ]

    params <- .fit_circle_on_3_points(points_subset)

    # Compute residuals for all points
    distances <- sqrt((points[, 1] - params[1])^2 + (points[, 2] - params[2])^2)
    residuals <- abs(distances - params[3])

    # Count inliers (points whose residuals are below the threshold)
    inliers <- sum(residuals < inlier_threshold)

    # Update best model if more inliers are found
    if (inliers > max_inliers)
    {
      max_inliers <- inliers
      best_circle <- params
    }
  }

  if (is.null(best_circle))
  {
    center_x = mean(points[,1])
    center_y = mean(points[,2])
    radius = 0.001
    z = mean(z)
    angle_range_degrees = 0
    inlier_ratio = 0
    percentage_inside = 0
    inliers = rep(TRUE, nrow(points))
    cfqi = 0
  }
  else
  {

    center_x <- best_circle[1]
    center_y <- best_circle[2]
    radius   <- best_circle[3]

    # Goodness of fit
    distances <- sqrt((points[, 1] - center_x)^2 + (points[, 2] - center_y)^2)
    residuals <- abs(distances - radius)
    inliers <- residuals < inlier_threshold

    # Angular range
    angle_res = 3
    angles <- atan2(points[inliers, 2] - center_y, points[inliers, 1] - center_x)
    angles <- ifelse(angles < 0, angles + 2 * pi, angles)
    angles <- sort(angles*180/pi)
    rangles <- unique(round(angles/angle_res)*angle_res)

    angle_range_degrees = sum(diff(rangles) <= angle_res)*angle_res
    arc_score <- angle_range_degrees / 360

    # Inlier ratio
    inlier_ratio <- sum(inliers) / nrow(points)

    # Percentage of point inside the circle
    ninside = sum(distances < radius - inlier_threshold)
    percentage_inside = ninside/nrow(points)
    score_inside = 1-percentage_inside

    # CFQI calculation
    w_arc <- 0.4
    w_inside = 0.1
    w_inlier <- 0.4
    cfqi <- (w_arc * arc_score + w_inlier * inlier_ratio + w_inside / score_inside) |> as.numeric()
  }

  return(list(center_x = center_x,
              center_y = center_y,
              radius = radius,
              z = mean(z),
              covered_arc_degree = angle_range_degrees,
              percentage_inlier = inlier_ratio,
              percentage_inside = percentage_inside,
              inliers = which(inliers),
              CFQI = cfqi))
}

is.valid.circle = function(radius, angle_range, pinliner)
{
  if (radius < 0.02)  return(TRUE)
  if (radius < 0.05)  return(angle_range > 180 & pinliner > 30)
  if (radius < 0.12)  return(angle_range > 90 & pinliner > 50)
  return(angle_range > 180 & pinliner > 30)
}


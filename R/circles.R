fit_circle_on_3_points <- function(points_subset)
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

ransac_circle <- function(points, num_iterations = 100, inlier_threshold = 0.01)
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

    params <- fit_circle_on_3_points(points_subset)

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

#' @importFrom methods is
is.valid.circle = function(radius, angle_range, pinliner)
{
  if (radius < 0.02)  return(TRUE)
  if (radius < 0.05)  return(angle_range > 180 & pinliner > 30)
  if (radius < 0.12)  return(angle_range > 90 & pinliner > 50)
  return(angle_range > 180 & pinliner > 30)
}

add_circles3d <- function(x, center_x, center_y, radius, height)
{
  theta <- seq(0, 2 * pi, length.out = 50)
  cos_theta <- cos(theta)
  sin_theta <- sin(theta)

  for (i in seq_along(center_x))
  {
    xc <- center_x[i]
    yc <- center_y[i]
    r <- radius[i]
    h <- height[i]

    xx <- xc - x[1] + r * cos_theta
    yy <- yc - x[2] + r * sin_theta
    zz <- rep(h, length(theta))

    rgl::lines3d(xx, yy, zz, col = "red", lwd = 2)
  }
}


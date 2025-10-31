#' @export
find_seeds <- function(las, params)
{
  foliage <- clusterID <- max_diameter <- passage <- hag <- NULL

  # The point cloud is supposed to have passage hag and foliage
  attributes <- names(las)
  stopifnot("hag" %in% attributes)
  stopifnot("passage" %in% attributes)
  stopifnot("foliage" %in% attributes)

  # Extract some slices of wood (thikcness 3cm)
  heights  <- c(min(las$hag), params$seed$slice_at)
  thick    <- params$seed$slice_thickness
  sor_k    <- params$seed$sor_k
  sor_m    <- params$seed$sor_m
  slices   <- slice_poi(las, heights, thick)
  somewood <- lidR::filter_poi(slices, foliage == 0)
  somewood <- lidR::classify_noise(somewood, lidR::sor(k = sor_k, m = sor_m)) # very very aggressive sor
  somewood <- lidR::remove_noise(somewood)
  somewood@data <- somewood@data[, .(X,Y,Z, passage, hag)]

  # Connect the point into clusters
  seed <- lidR::connected_components(somewood, 0.05, 10, connectivity = 26)
  seed <- lidR::filter_poi(seed, clusterID != 0)

  # For each cluster search for circles. If we have a nice circle we have a tree
  is.valid.circle <- function(radius, angle_range, pinliner, pinside)
  {
    if (radius < 0.02)  return(TRUE)
    if (radius  < 0.05)  return(angle_range > 180 & pinliner > 30)
    if (pinside > 20)   return(FALSE)
    if (radius < 0.10)  return(angle_range > 90 & pinliner > 50)
    return(angle_range > 180 & pinliner > 30)
  }
  fit_circle_to_seed <- function(id)
  {
    cl <- seed[seed$clusterID == id]
    if (lidR::npoints(cl) < 10) return(NULL)
    circle <- ransac_circle(cl, num_iterations = 400, inlier_threshold = 0.02)

    valid  <- is.valid.circle(circle$radius, circle$covered_arc_degree, circle$percentage_inlier*100, circle$percentage_inside*100)

    if (FALSE)
    {
      if (valid) col = "darkgreen" else col = "red"
      plot(cl$X, cl$Y, asp = 1, main = id)
      symbols(circle$center_x, circle$center_y, circles = circle$radius, inches = FALSE, add = TRUE, fg = col)
      symbols(circle$center_x, circle$center_y, circles = circle$radius+0.01, inches = FALSE, add = TRUE, fg = col)
      symbols(circle$center_x, circle$center_y, circles = circle$radius-0.01, inches = FALSE, add = TRUE, fg = col)
    }

    if (valid) return(data.frame(X = circle$center_x, Y = circle$center_y, Z = circle$z, R = circle$radius, id = id))
    else return(NULL)
  }
  circles <- lapply(unique(seed$clusterID), fit_circle_to_seed)
  circles <- do.call(rbind, circles)

  if (FALSE)
  {
    x <- plot(seed, color = 'clusterID', pal = pastel.colors(500))
    for (i in 1:nrow(circles))
      add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], circles$Z[i])
  }

  # For each circle we exclude wood in a safe zone beyond the circles.
  # This allow to clean false positive around important trees and prevent
  # dummy connection caused by noise in the next connected component step
  px <- somewood$X
  py <- somewood$Y
  pz <- somewood$Z
  rm <- rep(FALSE, lidR::npoints(somewood))
  safe_zone <- 0.2
  for (i in 1:nrow(circles))
  {
    cx <- circles$X[i]
    cy <- circles$Y[i]
    cz <- circles$Z[i]
    r  <- circles$R[i]
    d  <- sqrt((px-cx)^2 + (py-cy)^2 +(pz-cz)^2)
    rm[d > (r + 0.02) & d  < (r + safe_zone)] <- TRUE
  }
  somewood <- somewood[!rm]

  # Keep point with passage. And lower than max slicing
  # Increasing min passage reduce branches included in passage, But
  # we keep any passage close to the ground to segment bushes
  th <- max(heights) + 0.1
  min_passage     <- params$seed$min_passage
  tree_passages   <- lidR::filter_poi(las, passage > min_passage, hag < th)
  bushes_passages <- lidR::filter_poi(las, passage > 1, hag < min(hag) + 0.5, passage < 10)
  tree_passages@data <- tree_passages@data[, .(X,Y,Z, passage, hag)]
  bushes_passages@data <- bushes_passages@data[, .(X,Y,Z, passage, hag)]
  all_passages   <- suppressWarnings(rbind(tree_passages, bushes_passages))

  if (FALSE)
  {
    x = plot(tree_passages, pal = "green")
    plot(bushes_passages, add = x, pal = "red")
    plot(somewood, add = x, pal = foliage.colors[1])
    for (i in 1:nrow(circles))
      add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], circles$Z[i])
  }

  if (FALSE)
  {
    x = plot(all_passages, pal = "green")
    plot(somewood, add = x, pal = foliage.colors[1])
    for (i in 1:nrow(circles))
      add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], circles$Z[i])
  }

  # Convert to points
  generate_circle_points <- function(x, y, z, r, step = 0.1)
  {
    # Circumference
    circumference <- 2 * pi * r
    # Number of points for approximately every `step` meters
    n_points <- ceiling(circumference / step)
    # Angles for points
    theta <- seq(0, 2*pi, length.out = n_points + 1)[-1]  # remove last point to avoid duplicate
    # Generate points
    data.frame(
      X = x + r * cos(theta),
      Y = y + r * sin(theta),
      Z = rep(z, n_points)
    )
  }
  circle_points_list <- lapply(1:nrow(circles), function(i) {
    generate_circle_points(circles$X[i], circles$Y[i], circles$Z[i], circles$R[i], step = 0.05)
  })
  circle_points = do.call(rbind, circle_points_list)
  circle_points$passage = 1000
  circle_points$hag = 0
  circle_points = suppressWarnings(LAS(circle_points, header(las)))


  # Bind the wood and the passage in order to compute
  # connected component and merge passage from the same trees

  temp = suppressWarnings(rbind(somewood, tree_passages, bushes_passages, circle_points))
  res    <- round(params$decimation$barycentric_predecimation_resolution*1.5, 2)
  temp$Z <- temp$Z * 0.8
  temp   <- lidR::connected_components(temp, res, 1, name = "treeID", connectivity = 26)
  temp$Z <- temp$Z / 0.8

  if (FALSE)
  {
    x <- plot(temp, color = "treeID")
    #plot(passage, add = x)
    for (i in 1:nrow(circles))
      add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], circles$Z[i])
  }

  # Retain only the seed below BH
  seeds <- lidR::filter_poi(temp, passage > 0)
  seeds <- lidR::filter_poi(seeds, hag < 1)
  seeds
}

slice_poi = function(las, heights, thinkness = 0.02)
{
  # Build dynamic filter for slices
  slice_filter <- Reduce(`|`, lapply(heights, function(s)
  {
    (las$hag > (s-thinkness/2) & las$hag < (s + thinkness/2))
  }))

  lidR::filter_poi(las, slice_filter)
}

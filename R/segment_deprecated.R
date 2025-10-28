#' Extend the trees to the DTM
#'
#' Extend the trees to the DTM. It uses fit circle produced by \link{measure_diameters} and remove,
#' for each tree, every point below this circle and replace by a synthetic cylinder. This ensure a
#' proper cleaning of the bottom of the tree, the extension to the ground and a more robust QSM.
#'
#' @param trees LAS object from lidR
#' @param dtm SpatRaster from terra. The DTM.
#' @param circles data.frame produced by \link{measure_diameters}
#' @param extra_height The trees can be prolonged below the ground.
#' @noRd
tree_extensions = function(trees, dtm, circles, extra_height = 0.15)
{
  pointID <- treeID <- Z <- . <- NULL

  extensions = list()
  trees@data$pointID = 1:lidR::npoints(trees)

  data.table::setindex(trees@data, treeID)

  # Initialize progress bar
  pb <- utils::txtProgressBar(min = 0, max = nrow(circles), style = 3)

  for (i in 1:nrow(circles))
  {
    id  = circles$treeID[i]
    x = circles$center_x[i]
    y = circles$center_y[i]
    z = circles$center_z[i]
    r = circles$radius[i]
    HAG = circles$center_hag[i]
    circle = list(center_x = x, center_y = y, z = z, radius = r)

    loc = matrix(c(x, y), ncol = 2)
    zgnd = terra::extract(dtm, loc, method = "bilinear")
    hcyl = as.numeric(z - zgnd + extra_height)

    ids = trees@data[.(id), on = "treeID"][Z < z]$pointID
    trees$pointID[ids] = 0

    extension = generate_cylinder_points(circle, height = hcyl)
    lidR::quantize(extension$X, trees@header[["X scale factor"]], trees@header[["X offset"]])
    lidR::quantize(extension$Y, trees@header[["Y scale factor"]], trees@header[["Y offset"]])
    lidR::quantize(extension$Z, trees@header[["Z scale factor"]], trees@header[["Z offset"]])
    extension$treeID = id
    extensions[[as.character(id)]] = extension

    # Update progress bar
    utils::setTxtProgressBar(pb, i)
  }

  # Close progress bar
  close(pb)

  extensions = do.call(rbind, extensions)
  extensions$anisotropy = 1
  extensions$wood = TRUE
  data.table::setDT(extensions)

  extensions = lidR::LAS(extensions)

  trees = lidR::filter_poi(trees, pointID > 0)
  trees@data$pointID = NULL
  trees = weld_extension(trees, extensions)

  return(trees)
}


weld_extension <- function(trees, extensions, fill_value = 0L)
{
  ..all_cols <- NULL

  df1 = trees@data
  df2 = extensions@data

  all_cols <- union(names(df1), names(df2))

  # Add missing columns with the fill value
  for (col in setdiff(all_cols, names(df1))) {
    df1[[col]] <- fill_value
  }
  for (col in setdiff(all_cols, names(df2))) {
    df2[[col]] <- fill_value
  }

  # Ensure column order matches
  df1 <- df1[, ..all_cols]
  df2 <- df2[, ..all_cols]

  # Combine rows
  trees@data = rbind(df1, df2)
  trees = lidR::las_update(trees)

  return(trees)
}


#' Find a seed for each tree
#'
#' To segment individual tree instances, the first step is to identify seed points for each tree.
#' This is achieved by extracting a slice of points near the ground. The wood/foliage
#' semantic segmentation must be performed first using \link{segment_foliage} to ensure
#' that wood points are correctly identified. The method uses only points that are labelled 'wood'.
#' The method then applies connected component clustering followed by RANSAC circle fitting to reliably
#' assign a single seed to each tree.
#'
#' @param las A LAS object from lidR.
#' @param ... Unused. Additional parameters beyond \code{...} should generally remain unchanged,
#'   except in edge cases.
#' @param max_diameter Maximum expected tree diameter (in meters). Used to filter out invalid
#'   RANSAC-fitted circles.
#' @param slice_seeds_at A numeric vector of two values defining the height range (in meters)
#'   for slicing the point cloud.
#' @param res Resolution for connected component clustering.
# @param smooth Smoothing radius (in meters) applied to the slice to improve the clustering process.
#' @param tree_spacing numeric. If the scene is a plantation, providing the spacing of the trees
#' guarantee much better seeds finding by providing prior information. We now know that trees are
#' regularly spaced.
#' @noRd
find_seeds_old = function(las, slice_seeds_at = c(0.5, 0.8), ..., max_diameter = 0.6, res = 0.025, tree_spacing = NULL)
{
  treeID <- clusterID <- X <- Y <- Z <-  hag <- hag_max <- hag_min <- foliage <- NULL

  if (!is.null(tree_spacing) && tree_spacing > 0)
  {
    res = tree_spacing/4
  }

  # The point cloud must have hag, anisotropy and foliage computed
  attributes = names(las)
  stopifnot("foliage" %in% attributes)
  stopifnot("hag" %in% attributes)

  cat("Finding seeds with a clusetring approach\n")

  seed = lidR::filter_poi(las, hag > slice_seeds_at[1], hag < slice_seeds_at[2], foliage == FALSE)
  seed = smooth3d(seed, 0.03)
  seed$Z = seed$Z * 0.01
  seed = lidR::connected_components(seed, res, 10)
  seed = lidR::filter_poi(seed, clusterID != 0)
  seed$Z = seed$Z * 100
  #plot(seed, color = "clusterID", pal = pastel.colors(500)) |> add_dtm3d(dtm)

  fit_circle_to_seed = function(id, max_radius)
  {
    cl = seed[seed$clusterID == id]
    if (lidR::npoints(cl) < 10) return(NULL)
    circle = ransac_circle(cl, num_iterations = 400)
    valid  = is.valid.circle(circle$radius, circle$covered_arc_degree, circle$percentage_inlier*100, circle$percentage_inside*100)
    if (valid) return(data.frame(X = circle$center_x, Y = circle$center_y, Z = circle$z, R = circle$radius, id = id))
    else return(NULL)
  }

  cat("Fitting RANSAC circles to each seed\n")

  circles = lapply(unique(seed$clusterID), fit_circle_to_seed, max_radius = max_diameter/2)
  circles = do.call(rbind, circles)

  sfcircles = NULL
  if (!is.null(circles))
  {
    sfcenters = sf::st_as_sf(circles, coords = c("X", "Y", "Z"))
    sfcircles = sf::st_buffer(sfcenters, circles$R*1.20, nQuadSegs = 10)
  }

  f = function(x,y,z)
  {
    Z = c(round(min(z), 3), round(max(z), 3))
    dZ = diff(Z)

    bottom = z < Z[2] - dZ/2
    top = z > Z[1] + dZ/2

    X = c(round(stats::median(x[bottom]), 3), round(stats::median(x[top]), 3))
    Y = c(round(stats::median(y[bottom]), 3), round(stats::median(y[top]), 3))

    return(list(X = X, Y = Y,Z = Z))
  }

  cat("Seed correction with RANSAC circles\n")

  seeds = seed@data[, f(X,Y,Z), by = clusterID]
  seeds = stats::na.omit(seeds)

  sfseeds = sf::st_as_sf(seeds, coords = c("X", "Y", "Z"))

  extraseed = NULL

  if (!is.null(sfcircles))
  {
    intersect = sf::st_intersects(sfcircles, sfseeds)
    ii = lapply(intersect, length)
    ii = which(ii > 2)

    for (iii in ii)
    {
      ids = sfseeds[intersect[[iii]],]$clusterID
      sfseeds$clusterID[sfseeds$clusterID %in% ids] = ids[1]
    }

    cat("Seed generation from RANSAC circles\n")

    sfcenters = sfcenters[sfcenters$id %in% sfseeds$clusterID,]
    extraseed = sf::st_buffer(sfcenters, sfcenters$R*0.9, nQuadSegs = 3)
    extraseed$Z = sf::st_coordinates(sfcenters)[,3]
    extraseed = sf::st_cast(extraseed, "POINT")
    extraseed$R = NULL
    coord = as.data.frame(sf::st_coordinates(extraseed))
    coord$Z = extraseed$Z
    coord$clusterID = extraseed$id
    extraseed = sf::st_as_sf(coord, coords = c("X", "Y", "Z"))
  }

  fullseed = rbind(sfseeds, extraseed)



  #circ = circles[61,]
  #see = sfseeds[intersect[[61]],]
  #plot(sf::st_geometry(circ), axes = T)
  #plot(see, add = T, pch = as.numeric(as.factor(see$clusterID)), col = "black")
  #plot(seed[seed$clusterID %in% see$clusterID], color = "clusterID")

  names(fullseed)[1] = "treeID"

  #plot(seed, color = "clusterID", pal = pastel.colors(500)) |> add_dtm3d(dtm) |> add_treetops3d(sfseeds, radius = 0.05)

  return(fullseed)
}

#' Post production tool
#'
#' Apply this function to fix some issue in the instance segmentation. If you have no issue but apply
#' this function this should not have any effect. This function fits circles on each tree at different
#' levels and tries to figure out if two segmented trees are actually the same one badly segmented because
#' of bad seeds. In some cases tt is possible that one tree have two different seeds leading to two
#' different instance. This function tries to fix that.
#'
#' @param las LAS object from lidR
#' @param max_height maximum height above the lowest point of a tree to fit circles
#' @param slice_thickness slice thickness to fit circles
#' @param max_diameter a threshold to protect against bad circle fitting. If a circle is bigger than
#' that it is not valid. Set it the maximum expected diameter of a tree.
#' @noRd
fix_split_trees = function(las, max_height = 4, slice_thickness = 0.25, max_diameter = 0.6)
{
  treeID <- X <- Y <- Z <-  hag <- hag_max <- hag_min <- foliage <- NULL

  attributes = names(las)
  stopifnot("foliage" %in% attributes)
  stopifnot("hag" %in% attributes)

  max_radius = max_diameter/2

  # Define a function to fit circle
  fit_circle_to_seed = function(id)
  {
    #cat("Tree", id, ":")
    keep = slice$treeID == id
    keep[is.na(keep)] = FALSE
    cl = slice[keep]
    if (lidR::npoints(cl) < 10)
    {
      #cat(" not enought points\n")
      return(NULL)
    }
    circle = ransac_circle(cl, num_iterations = 400, inlier_threshold = 0.02)
    valid = is.valid.circle(circle$radius, circle$covered_arc_degree, circle$percentage_inlier*100, circle$percentage_inside*100)
    if (valid)
    {
      #cat("\n")
      return(data.frame(X = circle$center_x, Y = circle$center_y, Z = circle$z, R = circle$radius, id = id))
    }
    else
    {
      #cat(" invalid\n")
      return(NULL)
    }
  }

  zmin = min(las$hag)
  offsets = seq(0, max_height,slice_thickness)

  # For each slice apply the fit to each tree
  # if two circle intersect they are from the same tree
  for (k in 1:length(offsets))
  {
    slice = lidR::filter_poi(las, hag > zmin+offsets[k], hag < zmin+offsets[k]+slice_thickness, foliage == FALSE)

    #x = plot(slice, color = "treeID")
    circles = lapply(unique(slice$treeID), fit_circle_to_seed)
    circles = do.call(rbind, circles)

    if (is.null(circles)) next

    #for (i in 1:nrow(circles)) add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], 0)

    circles = sf::st_as_sf(circles, coords = c("X", "Y"))
    circles = sf::st_buffer(circles, circles$R*1.20)

    intersect = sf::st_intersects(circles)
    intersect = Filter(function(x) length(x) > 1, intersect)
    intersect = unique(intersect)

    if (length(intersect) == 0) next

    for (i in intersect)
    {
      intersecting_circles = circles[i,]
      ids = intersecting_circles$id

      intersecting_circles = sf::st_geometry(intersecting_circles)

      if (length(intersecting_circles) > 2)
        warning("More than 2 rings intersecting. This case is not handled yet.")

      area1 = mean(sf::st_area(intersecting_circles))
      area2 = sf::st_area(sf::st_intersection(intersecting_circles[1],intersecting_circles[2]))


      if (area2 > area1*0.5)
      {
        #plot(filter_poi(las, treeID %in% ids), color = "treeID", axis = T)
        #plot(filter_poi(las, treeID %in% ids), color = "foliage")
        cat("Combining trees", ids, "into tree", ids[1], "\n")
        las$treeID[las$treeID %in% ids] = ids[1]
      }
    }
  }

  return(las)
}

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
#' @noRd
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
    tt = lidR::decimate_points(tt, lidR::random_per_voxel(0.02))


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

generate_cylinder_points <- function(circle, height = 0.5, n_points = 15000)
{
  # Extract circle parameters
  center_x <- circle$center_x
  center_y <- circle$center_y
  radius <- circle$radius
  z_top <- circle$z
  z_bottom <- z_top - height

  # Create a grid of points
  cylinder_points <- data.frame(
    X = 0,
    Y = 0,
    theta = stats::runif(n_points, 0, 2*pi),
    Z = stats::runif(n_points, z_bottom, z_top)
  )

  # Convert to Cartesian coordinates
  cylinder_points$X <- center_x + radius * cos(cylinder_points$theta)
  cylinder_points$Y <- center_y + radius * sin(cylinder_points$theta)

  # Drop theta column (optional)
  cylinder_points$theta <- NULL

  # Return the result as a data frame
  return(cylinder_points)
}

align_to_z <- function(main_axis)
{
  # Normalize the main axis vector
  main_axis <- main_axis / sqrt(sum(main_axis^2))

  # Z-axis unit vector
  z_axis <- c(0, 0, 1)

  # Cross product to find the rotation axis
  rotation_axis <- c(
    main_axis[2] * z_axis[3] - main_axis[3] * z_axis[2],
    main_axis[3] * z_axis[1] - main_axis[1] * z_axis[3],
    main_axis[1] * z_axis[2] - main_axis[2] * z_axis[1]
  )

  # Normalize the rotation axis
  axis_length <- sqrt(sum(rotation_axis^2))
  if (axis_length > 1e-6) { # Avoid division by zero
    rotation_axis <- rotation_axis / axis_length
  } else {
    # If the main axis is already aligned with Z, return identity matrix
    return(diag(3))
  }

  # Compute the angle between main_axis and Z-axis
  angle <- acos(sum(main_axis * z_axis))

  # Construct the rotation matrix using Rodrigues' rotation formula
  K <- matrix(c(
    0, -rotation_axis[3], rotation_axis[2],
    rotation_axis[3], 0, -rotation_axis[1],
    -rotation_axis[2], rotation_axis[1], 0
  ), nrow = 3, byrow = TRUE)

  R <- diag(3) + sin(angle) * K + (1 - cos(angle)) * (K %*% K)

  return(R)
}

combine_with_fill <- function(df1, df2, fill_value = 0L)
{
  ..all_cols <- NULL

  all_cols <- union(names(df1), names(df2))

  # Add missing columns with the fill value
  for (col in setdiff(all_cols, names(df1))) {
    df1[[col]] <- fill_value
  }
  for (col in setdiff(all_cols, names(df2))) {
    df2[[col]] <- fill_value
  }

  # Ensure column order matches
  df1 <- df1[, ..all_cols]
  df2 <- df2[, ..all_cols]

  # Combine rows
  rbind(df1, df2)
}

smooth3d = function(las, radius = 0.05, weight = NULL, ncpu = 8, progress = TRUE)
{
  ti = Sys.time()

  if (is.function(radius)) radius = radius(las$Z)
  if (length(radius) == 1L) radius = rep(radius, lidR::npoints(las))
  stopifnot(length(radius) == lidR::npoints(las))
  stopifnot(all(radius > 0))
  if (is.null(weight)) weight = rep(1, lidR::npoints(las))
  stopifnot(length(weight) == lidR::npoints(las), radius > 0)

  z = cpp_smooth3d(las, radius = radius, weight = weight, ncpu = ncpu, progress, FALSE)

  tf = Sys.time()
  print(difftime(tf,ti))

  las$X = z$X
  las$Y = z$Y
  las$Z = z$Z
  las
}
compute_point_network <- function(dec, k, max_gap = 1, wood_mask = NULL, cost_factors = NULL, power = 3, downward = FALSE)
{
  if (!is.null(wood_mask))
  {
    stopifnot(
      !is.null(cost_factors))
  }

  if (!is.null(cost_factors))
  {
    stopifnot(
      !is.null(cost_factors$wood2wood),
      !is.null(cost_factors$wood2leaf),
      !is.null(cost_factors$leaf2leaf))
  }

  net <- compute_network(dec, k = k)

  # Gaps have an infinite cost
  net$cost[net$cost > max_gap] <- Inf

  # It is more expensive to move in large steps
  net$cost <- net$cost^power

  # It is more expensive to move downward. The path finder starts from
  # a seed below the ground. If should reach any target points by moving upward
  # preferentially. If it moves downward it is not following a tree. Maybe a
  # branch bending. It is ok but expensive.
  X1 <- dec$X[net$from]; X2 <- dec$X[net$to]
  Y1 <- dec$Y[net$from]; Y2 <- dec$Y[net$to]
  Z1 <- dec$Z[net$from]; Z2 <- dec$Z[net$to]
  dz <- Z1-Z2
  magnitude <- sqrt((X1-X2)^2 + (Y1-Y2)^2 + dz^2)
  cos_theta <- -dz / magnitude

  if (downward) cos_theta = -cos_theta

  angle_degree <- acos(pmin(pmax(cos_theta, -1), 1)) * 180/pi

  f <- function(x){ y = exp(log(100)/100*x); y[x>100]=100; y }
  net$cost <- net$cost * f(angle_degree)

  # Optional wood/leaf adjustment costs
  if (!is.null(wood_mask) & !is.null(cost_factors))
  {
    is_wood1 <- !wood_mask[net$from]
    is_wood2 <- !wood_mask[net$to]
    wood2wood <- is_wood1 & is_wood2
    leaf2leaf <- !is_wood1 & !is_wood2
    wood2leaf <- is_wood1 & !is_wood2
    net$cost[wood2wood] = net$cost[wood2wood] * cost_factors$wood2wood
    net$cost[leaf2leaf] = net$cost[leaf2leaf] * cost_factors$leaf2leaf
    net$cost[wood2leaf] = net$cost[wood2leaf] * cost_factors$wood2leaf
  }

  return(net)
}






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
#' @export
fix_split_trees = function(las, max_height = 2, slice_thickness = 0.25, max_diameter = 0.6)
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
    valid = is.valid.circle(circle$radius, circle$covered_arc_degree, circle$percentage_inlier*100)
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

    #for (i in 1:nrow(circles)) add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], 0)

    circles = sf::st_as_sf(circles, coords = c("X", "Y"))
    circles = sf::st_buffer(circles, circles$R*1.20)

    intersect = sf::st_intersects(circles)
    intersect = Filter(function(x) length(x) > 1, intersect)
    intersect = unique(intersect)

    if (length(intersect) == 0) next

    for (i in intersect)
    {
      print(i)
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

#' Post production tool
#'
#' Apply this function to fix some issue in the instance segmentation. If you have no issue but apply
#' this function this should not have any effect. This function look at the wood part of each trees
#' an test with a connected component algorithm if there are more than one big cluster. This can
#' happen when a seed is missing with \link{find_seeds}. In this case a few wood point may be assigned
#' to the wrong tree because this was the best seed to reach. Those point are reasigned to foliage class
#' in order to do not mess-up QSM software
#'
#' @param las LAS object from lidR with semantic and instance segmentation
#' @export
fix_small_isolated_low_clusters = function(las)
{
  treeID <- clusterID <- foliage <- hag <- foliage <- NULL

  las@data$pointID = 1:lidR::npoints(las)
  sub = lidR::filter_poi(las, foliage == FALSE, hag < 3)

  for (id in unique(las$treeID))
  {
    cat("Tree", id)
    tt = lidR::filter_poi(sub, treeID == id)
    if (lidR::is.empty(tt))
    {
      cat("\n")
      next
    }
    tt$Z = tt$Z * 0.1
    tt = lidR::connected_components(tt, 0.05, 200)
    tt$Z = tt$Z * 10

    ids = 1
    n = length(unique(tt$clusterID))
    if (n > 1)
    {
      ids = table(tt$clusterID)
      ids = as.numeric(names(ids[which.max(ids)]))
      cat(" :", n, "clusters found. Smaller cluster(s) re-assigned as foliage")
      pid = tt$pointID[tt$clusterID != ids]
      las$foliage[pid] = TRUE
      tt = lidR::filter_poi(tt, clusterID == ids)
      #plot(tt, color = "clusterID")
      #cat("  ", id, "\n")
      #plot(tt, color = "foliage", pal = c("chocolate4", "darkgreen"))
    }
    else
    {
      cat(" ok")
    }

    cat("\n")
  }

  las@data$pointID = NULL

  return(las)
}

generate_cylinder_points <- function(circle, height = 0.5, n_points = 2000)
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


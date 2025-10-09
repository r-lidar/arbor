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


#' Remove Small Trees and Clean Low Understory
#'
#' This function removes small trees that are likely to be poorly segmented and filters out
#' low understory vegetation based on height.
#'
#' @param las A LAS object from lidR.
#' @param max_height Trees with a height less than this threshold (in meters) will be removed.
#' @export
clean_small_cluster = function(las, max_height = 3)
{
  treeID <- hag <- hag_max <- hag_min <- NULL

  attributes = names(las)
  stopifnot("treeID" %in% attributes)
  stopifnot("hag" %in% attributes)

  ans = las@data[!is.na(treeID), list(hag_max = max(hag), hag_min = min(hag)), by = treeID]
  ans = ans[hag_max > max_height & hag_min < max_height]

  trees = lidR::filter_poi(las, treeID %in% ans$treeID & !is.na(treeID))
  trees
}


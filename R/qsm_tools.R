#' Calculate Distance from Point to Cylinder Surface
#'
#' Computes the distance from a point to the surface of a cylinder segment.
#'
#' @param point A numeric vector representing the point in 3D space.
#' @param start A numeric vector representing the start of the cylinder segment.
#' @param end A numeric vector representing the end of the cylinder segment.
#' @param radius A numeric value for the cylinder radius.
#' @return The shortest distance from the point to the cylinder surface.
#' @noRd
#' @md
point_to_cylinder_distance <- function(point, start, end, radius)
{
  infinity = rep(.Machine$double.xmax)

  # Cylinder's axis
  axis_vector <- end - start
  axis_length <- sqrt(sum(axis_vector^2))
  axis_direction <- axis_vector / axis_length

  # Project point onto the cylinder's axis
  point_vector <- point - start
  projection_length <- sum(point_vector * axis_direction)

  # Check if the projection is beyond the cylinder segment
  if (projection_length < 0)
    closest_point_on_axis <- infinity
  else if (projection_length > axis_length)
    closest_point_on_axis <- infinity
  else
    closest_point_on_axis <- start + projection_length * axis_direction

  # Calculate the distance from the point to the closest point on the axis
  distance_to_axis <- sqrt(sum((point - closest_point_on_axis)^2))
  return(max(0, distance_to_axis - radius))
}

#' Query Points Near a Cylinder
#'
#' Filters points from a LAS object that are within a cylinder's radius.
#'
#' @param las A LAS object representing point cloud data.
#' @param cyl A list containing cylinder parameters.
#' @return A filtered LAS object with points within the specified cylinder.
#' @noRd
#' @md
query_cylinder = function(las, cyl)
{
  tx = cyl$translateX
  ty = cyl$translateY
  tz = cyl$translatez
  start <- c(cyl$startX - tx, cyl$startY - ty, cyl$startZ - tz)
  end <- c(cyl$endX - tx, cyl$endY - ty, cyl$endZ- tz)
  r = cyl$radius

  u = sf::st_coordinates(las)

  d = sapply(1:nrow(u), function(i)
  {
    p = u[i,]
    d = point_to_cylinder_distance(p, start, end, r)
    d
  })

  sub = lidR::filter_poi(las, d <= r+0.1)
  return(sub)
}

#' Calculate Volume Metrics for QSM Data
#'
#' Summarizes volume metrics in QSM data, plotting cumulative volumes and attributes.
#'
#' @param qsm A data frame of QSM data.
#' @param ... passed to qsm_outlier_detection()
#' @param plot boolean display graphics
#' @return The total volume of the QSM object.
#' @export
#' @md
volume = function(qsm, ..., plot = F)
{
  stopifnot("volume" %in% names(qsm))
  V = sum(qsm$volume)
  return(V)
}

#' Calculate Distance to Root for QSM Segments
#'
#' Computes the distance from each segment to the root of the QSM structure.
#'
#' @param qsm A data frame of QSM segment data.
#' @return The modified QSM data frame with calculated distances to root.
#' @noRd
#' @md
qsm_distance_to_root = function(qsm)
{
  qsm$dist2root <- NA_real_
  qsm$dist2root[qsm$parent_ID == 0] <- 0  # Set root distance to 0

  cylinder_id_queue <- list(1)  # start from the root (assuming ID 1)

  # Distance to root
  while (length(cylinder_id_queue) > 0)
  {
    currentID <- cylinder_id_queue[[1]]
    currentIndex <- which(qsm$cyl_ID == currentID)

    parentID <- qsm$parent_ID[currentIndex]
    if (parentID != 0)
    {
      parentIndex <- which(qsm$cyl_ID == parentID)
      qsm$dist2root[currentIndex] <- qsm$dist2root[parentIndex] + qsm$length[currentIndex]
    }

    # Enqueue all children of the current node
    childIDs <- qsm$cyl_ID[qsm$parent_ID == currentID]
    cylinder_id_queue <- c(cylinder_id_queue[-1], as.list(childIDs))
  }

  return(qsm)
}

filter_tree = function(tree)
{
  attributes = names(tree)

  if ("foliage" %in% attributes)
  {
    tree = lidR::filter_poi(tree, foliage == FALSE)
  }

  if ("treeID" %in% attributes)
  {
    if (length(unique(tree$treeID)) != 1)
      stop("The point cloud must contain a single tree", call. = FALSE)
  }

  return(tree)
}

compute_rotation_matrix <- function(start, end)
{
  # Vector representing the axis
  v <- end - start
  v <- v / norm(v, type = "2")  # Normalize

  # Target vector (Z-axis)
  z <- c(0, 0, 1)

  # Check if already aligned
  if (all(abs(v - z) < 1e-8)) return(diag(3))
  if (all(abs(v + z) < 1e-8)) return(-diag(3))

  # Rotation axis and angle
  axis <- pracma::cross(v, z)
  axis <- axis / norm(axis, type = "2")
  angle <- acos(pracma::dot(v, z))

  # Rodrigues' rotation formula
  K <- matrix(c(0, -axis[3], axis[2], axis[3], 0, -axis[1], -axis[2], axis[1], 0), 3, 3, byrow = TRUE)
  R <- diag(3) + sin(angle) * K + (1 - cos(angle)) * (K %*% K)
  return(R)
}





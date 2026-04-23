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
  qsm$length = with(qsm, sqrt((startX - endX)^2 +  (startY - endY)^2 +  (startZ - endZ)^2))
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
  foliage <- NULL
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





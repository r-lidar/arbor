#' Simple tools for QSM
#'
#' Simple tools for QSM.
#'
#' @param qsm A QSM
#' @export
#' @rdname qsm_tools
#' @md
qsm_volume = function(qsm)
{
  stopifnot("volume" %in% names(qsm))
  V = sum(qsm$volume)
  return(V)
}

#' @rdname qsm_tools
#' @export
qsm_message = function(qsm)
{
  msg = attr(qsm, "message")
  msg
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

#' Add Synthetic Ground Points to a Single-Tree Point Cloud
#'
#' Adds a set of synthetic ground points below a single-tree point cloud. Useful to perform
#' arbor operations on already externally isolated trees
#'
#' @param las A `LAS` object containing a single-tree point cloud.
#' @param n integer. Number of ground points
#'
#' @return A `LAS` object containing the original points plus n synthetic
#' ground points.
#'
#' @export
add_single_tree_ground = function(las, n = 1000)
{
  z   <- min(las$Z)
  bb  <- lidR::st_bbox(las)
  xg  <- stats::runif(n, bb[1]-1, bb[3]+1)
  yg  <- stats::runif(n, bb[2]-1, bb[4]+1)
  lidR::quantize(xg, 0.001, las@header[["X offset"]])
  lidR::quantize(xg, 0.001, las@header[["Y offset"]])
  gnd <- data.frame(X = xg, Y = yg, Z = z)
  gnd <- suppressWarnings(lidR::LAS(gnd, header = las@header))
  suppressWarnings(rbind(las, gnd))
}






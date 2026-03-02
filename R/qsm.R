#' Generate a QSM from a single tree point cloud
#'
#' This function processes a tree point cloud to generate a Quantitative Structure Model (QSM) using
#' the bottom to top approach to build the skeleton. Then if fits a radius to each edge of the skeleton
#' with a ransac circle fitting. Last it applies a allometric model to correct irrelevant cylinder
#' and compute a more organic tree. See the [Arbor book] for mode details.
#'
#' @param tree A `LAS` object from lidR containing a single tree point cloud. Only
#' the point labelled as wood will be used for QSM. Can alternatively be a string to a file
#' @param step The slice thinkness during the slicing operation (keep as is)
#' @param cl_dist Passed to eps in dbscan::dbscan (keep as is)
#' @param max_d Maximum distance between cluster (keep as is)
#' @param apex last radius of the branch tips
#' @param ... unused
#' @examples
#' f <- system.file("extdata", "tree_qsm.laz", package="arbor")
#' tree <- lidR::readLAS(f)
#' qsm = qsm(tree)
#' x = plot_semantic(tree)
#' plot_qsm(qsm, add = x, color = "branch_order", cylinder = TRUE)
#' @export
#' @seealso \link{qsm_write} \link{qsm_read} \link{qsm_dbh} \link{qsm_stats}
qsm = function(tree, step = 0.2, cl_dist = 0.1, max_d = 0.1, apex = 0.0025, ...)
{
  #step = 0.2; cl_dist = 0.1; max_d = 0.1; apex = 0.0025

  t0 <- tic() ; gc()

  if (is.character(tree))
  {
    if (!file.exists(tree))
      stop(paste0(tree, " does not exist"))

    tree = lidR::readLAS(tree, select = "0")
  }

  j    <- which.min(tree$Z)
  tx   <- tree$X[j]
  ty   <- tree$Y[j]
  tz   <- tree$Z[j]

  tree <- shift(tree, tx, ty, tz)   # Move to origin for numerical stability
  tree <- filter_tree(tree)         # Remove foliage
  tree <- clean_tree_butt(tree)     # R

  qsm  <- qsm_skeleton(tree, step, cl_dist, max_d) # c++
  qsm  <- qsm_architecture(qsm)                    # c++
  qsm  <- qsm_smooth(qsm, niter = 1)               # c++
  qsm  <- qsm_detect_weird_butt(qsm)               # R

  d    <- estimate_prolongation(tree, qsm)         # R

  qsm  <- qsm_prolongation(qsm, d)                 # R
  qsm  <- qsm_radius(qsm, tree, tip_radius = apex) # R & c++
  qsm  <- qsm_volume(qsm)
  qsm  <- shift(qsm, -tx, -ty, -tz)

  order = c("startX", "startY", "startZ", "endX", "endY", "endZ", "cyl_ID", "parent_ID", "axis_ID", "branch_order","subtree_length", "radius", "volume")
  data.table::setcolorder(qsm, order)

  qsm <- set_qsm_class(qsm)

  st_crs(qsm) <- sf::st_crs(tree)

  toc(t0, space = "")

  return(qsm)
}

clean_tree_butt = function(tree)
{
  logger("Cleaning tree butt")

  tree@data$pointID <- 1:lidR::npoints(tree)
  bottom <- tree[tree$Z < min(tree$Z) + 0.25]
  bottom <- connected_components(bottom, 0.05, 10, connectivity = 26)

  if (length(unique(bottom$clusterID)) > 1)
  {
    logger("Multiple clusters at the bottom of the tree detected. Automatic cleaning triggered.", level = "WARN")

    t <- table(bottom$clusterID)
    i <- as.numeric(names(which.max(t)))
    r <- bottom$pointID[bottom$clusterID != i]
    tree <- tree[-r]
  }

  return(tree)
}

estimate_prolongation = function(tree, qsm)
{
  # compute the prolongation length d
  d = 0
  if ("hag" %in% names(tree))
  {
    minZ <- min(qsm$startZ)
    hag  <- tree$hag[tree$Z <= minZ]
    d    <- max(hag)
  }

  return(d)
}

set_qsm_class <- function(x)
{
  stopifnot(data.table::is.data.table(x))
  class(x) <- c("qsm", class(x))
  x
}


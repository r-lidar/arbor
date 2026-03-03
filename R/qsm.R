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

  # Move to origin for numerical stability
  j    <- which.min(tree$Z)
  tx   <- tree$X[j]
  ty   <- tree$Y[j]
  tz   <- tree$Z[j]
  tree <- shift(tree, tx, ty, tz)

  tree <- filter_tree(tree)                        # R
  tree <- clean_tree_butt(tree)                    # c++

  qsm  <- qsm_skeleton(tree, step, cl_dist, max_d) # c++
  qsm  <- qsm_architecture(qsm)                    # c++
  qsm  <- qsm_smooth(qsm, niter = 1)               # c++
  qsm  <- qsm_detect_weird_butt(qsm)               # c++

  d    <- estimate_prolongation(tree, qsm)         # c++

  qsm  <- qsm_prolongation(qsm, d)                 # c++
  qsm  <- qsm_radius(qsm, tree, tip_radius = apex) # c++
  qsm  <- qsm_volume(qsm)
  qsm  <- shift(qsm, -tx, -ty, -tz)

  order = c("startX", "startY", "startZ", "endX", "endY", "endZ", "cyl_ID", "parent_ID", "axis_ID", "branch_order","subtree_length", "radius", "volume")
  data.table::setcolorder(qsm, order)

  qsm <- set_qsm_class(qsm)

  st_crs(qsm) <- sf::st_crs(tree)

  toc(t0, space = "")

  return(qsm)
}

set_qsm_class <- function(x)
{
  stopifnot(data.table::is.data.table(x))
  class(x) <- c("qsm", class(x))
  x
}


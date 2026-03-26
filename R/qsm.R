#' Generate a QSM from a single tree point cloud
#'
#' This function processes a tree point cloud to generate a Quantitative Structure Model (QSM) using
#' the bottom to top approach to build the skeleton. Then if fits a radius to each edge of the skeleton
#' with a ransac circle fitting. Last it applies a allometric model to correct irrelevant cylinder
#' and compute a more organic tree. See the [Arbor book] for mode details.
#'
#' @param tree A `LAS` object from lidR containing a single tree point cloud. Only
#' the point labelled as wood will be used for QSM. Can alternatively be a string to a file
#' @param params list See \link{parameters}.
#' @examples
#' f <- system.file("extdata", "tree_qsm.laz", package="arbor")
#' tree <- lidR::readLAS(f)
#' qsm = qsm(tree)
#' \dontrun{
#' x = plot_semantic(tree)
#' plot_qsm(qsm, add = x, color = "branch_order", cylinder = TRUE)
#' }
#' @export
#' @seealso \link{qsm_write} \link{qsm_read} \link{qsm_dbh} \link{qsm_stats}
qsm =  function(tree, params = default_arbor_parameters)
{
  params <- evaluate_penalty(params)
  qsm <- qsm_cpp(tree@data, params)
  qsm_finalize(qsm)
}

qsm_finalize = function(qsm)
{
  qsm <- qsm_volume(qsm)
  data.table::setDT(qsm)
  order <- c("startX", "startY", "startZ", "endX", "endY", "endZ", "cyl_ID", "parent_ID", "axis_ID", "branch_order","subtree_length", "radius", "volume")
  data.table::setcolorder(qsm, order)
  qsm <- set_qsm_class(qsm)
  qsm
}

set_qsm_class <- function(x)
{
  stopifnot(data.table::is.data.table(x))
  class(x) <- c("qsm", class(x))
  x
}


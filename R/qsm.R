#' Generate a QSM from a single tree point cloud
#'
#' This function processes a tree point cloud to generate a Quantitative Structure Model (QSM)
#' See the [Arbor book](<placeholder>) for mode details.
#'
#' @param tree A `LAS` object from lidR containing a single tree point cloud. Only
#' the point labelled as wood will be used for QSM.
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
#' @md
#' @seealso \link{qsm_write} \link{qsm_read} \link{qsm_dbh} \link{qsm_stats}
qsm =  function(tree, params = arbor_parameters_default)
{
  qsm <- qsm_cpp(tree@data, params)
  qsm_finalize(qsm)
}

qsm_finalize = function(qsm)
{
  data.table::setDT(qsm)
  order <- c("startX", "startY", "startZ", "endX", "endY", "endZ", "cyl_ID", "parent_ID", "axis_ID", "branch_order", "dist_to_root", "subtree_length", "radius")
  data.table::setcolorder(qsm, order)
  qsm <- as_qsm(qsm)
  msg = qsm_message(qsm)
  if (length(msg) > 0) warning(msg, call. = FALSE)
  qsm
}

as_qsm <- function(x)
{
  data.table::setDT(x)
  class(x) <- c("qsm", class(x))
  x
}


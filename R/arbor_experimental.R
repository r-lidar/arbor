#' Extract the local tree context from a LAS object
#'
#' Extracts points belonging to trees surrounding a target tree, based on
#' different spatial context definitions (contact, crown extent, or root zone).
#' This is useful for analysing neighbourhood interactions, competition, or
#' structural context around individual trees.
#'
#' @param las A \code{LAS} object (from \pkg{lidR}) containing a segmented forest
#'   point cloud with a \code{treeID} attribute.
#' @param tree Either a numeric tree ID corresponding to \code{treeID} in
#'   \code{las}, or a \code{LAS} object containing points from a single tree.
#' @param context Character string defining how the spatial context is computed.
#'   One of:
#'   \describe{
#'     \item{\code{"contact"}}{Keeps only neighbouring trees that are in direct
#'     spatial contact with the target tree, based on nearest-neighbour queries.}
#'     \item{\code{"extent"}}{Keeps all trees whose points intersect the convex
#'     hull of the target tree.}
#'     \item{\code{"root"}}{Keeps trees intersecting the convex hull of the target
#'     tree restricted to points below 2 m height above ground.}
#'   }
#' @param exclude_tree Logical. If \code{TRUE}, the target tree itself is removed
#'   from the returned context.
#'
#' @return A \code{LAS} object containing points from neighbouring trees that form
#'   the selected context.
#'
#' @examples
#' \dontrun{
#' las <- readLAS("segmented_forest.las")
#'
#' tree <- lidR::filter_poi(las, treeID == 10)
#' context <- extract_tree_context(las, tree, exclude_tree = TRUE)
#' x <- plot_instance(context)
#' plot(tree, add = x, pal = "red")
#' }
#'
#' @export
#' @family experimental
extract_tree_context = function(las, tree, exclude_tree = FALSE, verbose = TRUE)
{
  warn_experimental()

  # ---- Argument checks --------------------------------------------------

  if (!inherits(las, "LAS")) stop("'las' must be a LAS object")
  if (lidR::is.empty(las)) stop("'las' is empty")
  if (!"treeID" %in% names(las)) stop("'las' must contain a 'treeID' attribute")
  if (!is.logical(exclude_tree) || length(exclude_tree) != 1) stop("'exclude_tree' must be a single logical value")
  if (is.numeric(tree))
  {
    if (length(tree) != 1 || is.na(tree))
      stop("'tree' must be a single, non-missing numeric treeID")
  }

  if (is.numeric(tree))
  {
    id   <- tree
    tree <- lidR::filter_poi(las, treeID == id)
    if (is.empty(tree)) stop(paste0("No tree with treeID = ", id))
  }
  else
  {
    id   <- tree$treeID[1]
  }

  bb  <- sf::st_bbox(tree)
  roi <- lidR::clip_rectangle(las, bb[1], bb[2], bb[3], bb[4])
  ids <- unique(roi$treeID)
  ids <- na.omit(ids)

  if (exclude_tree)
  {
    ids <- ids[ids != id]
  }

  res <- lidR::filter_poi(las, treeID %in% ids)

  ll  <- lidR::decimate_points(res, lidR::random_per_voxel(0.25, 2))
  nn  <- lidR::knnx(ll, tree)
  idx <- ll$treeID[unique(as.integer(nn$nn.index))]
  idx <- unique(idx)
  res <- lidR::filter_poi(res, treeID %in% idx)

  return(res)
}


#' Semantic segmentation of a tree point cloud using a QSM
#'
#' Performs an experimental semantic segmentation of a single-tree point cloud
#' using distances derived from a Quantitative Structure Model (QSM). Points
#' are labelled as foliage or non-foliage based on their proximity to QSM
#' elements.
#'
#' @param las A \code{LAS} object containing points from a single tree
#' @param qsf A Quantitative Structure Forest associated with \code{las}.
#'
#' @family experimental
#' @export
qsf_segment_semantic = function(las, qsf)
{
  warn_experimental()
  if (inherits(qsf, "qsf")) qsf = data.table::rbindlist(qsf)
  res <- qsm_distances_cpp(qsf, las@data)
  b <- res$dist < res$radius*1.3 + 0.02
  las@data$foliage <- 1
  las@data$foliage[b] <- 0
  return(las)
}

warn_experimental = function()
{
  message("This is an experimental function that may not be retained in later versions of arbor.")
}

#' Generate a QSM from a single tree point cloud
#'
#' This function processes a tree point cloud to generate a Quantitative Structure Model (QSM) using
#' the bottom to top approach to build the skeleton. Then if fits a radius to each edge of the skeleton
#' either with a ransac circle fitting or a least square approach depending on if ransac successfully
#' fit a good circle. Last it applies an AdTree-like allometric model to correct irrelevant cylinder
#' and compute a more organic tree.
#'
#' @param tree A `LAS` object containing a single tree point cloud with required attributes. Only
#' the point labelled as wood will be used for QSM.
#' @param apex last radius of the branch tips
#'
#' @export
qsm = function(tree, step = 0.2, cl_dist = 0.1, max_d = 0.1, apex = 0.0025, ...)
{
  t0 = tic()
  #step = 0.2; cl_dist = 0.1; max_d = 0.1; apex = 0.005; power = 1.1; pure_model = FALSE; verbose = FALSE

  tree <- filter_tree(tree)

  # Move to origin for numerical stability
  j           <- which.min(tree$Z)
  tx          <- tree$X[j]
  ty          <- tree$Y[j]
  tz          <- tree$Z[j]
  tree@data$X <- tree@data$X-tx
  tree@data$Y <- tree@data$Y-ty
  tree@data$Z <- tree@data$Z-tz

  tree <- qsm_clean_tree_butt(tree)

  qsm <- qsm_skeleton(tree, step, cl_dist, max_d, verbose)
  qsm <- qsm_topology(qsm)
  qsm <- qsm_architecture(qsm)
  qsm <- qsm_detect_weird_butt(qsm)

  # compute the prolongation length d
  d = 0
  if ("hag" %in% names(tree))
  {
    minZ <- min(qsm$startZ)
    hag  <- tree$hag[tree$Z <= minZ]
    d    <- max(hag)
  }

  qsm <- qsm_prolongation(qsm, d)
  qsm <- qsm_radius(qsm, tree, R0, tip_radius = apex)
  qsm <- qsm_volume(qsm)

  qsm$startX  <- qsm$startX+tx
  qsm$startY  <- qsm$startY+ty
  qsm$startZ  <- qsm$startZ+tz
  qsm$endX    <- qsm$endX+tx
  qsm$endY    <- qsm$endY+ty
  qsm$endZ    <- qsm$endZ+tz

  toc(t0, space = "")

  return(qsm)
}

qsm_clean_tree_butt = function(tree)
{
  cat("Cleaning tree butt...") ; ti = tic()

  tree@data$pointID <- 1:lidR::npoints(tree)
  bottom <- tree[tree$Z < 0.25]
  bottom <- lidR::connected_components(bottom, 0.05, 10, connectivity = 26)

  if (length(unique(bottom$clusterID)) > 1)
  {
    cat("\n\033[33m  Multiple clusters at the bottom of the tree detected. Automatic cleaning triggered.\033[0m\n")

    t <- table(bottom$clusterID)
    i <- as.numeric(names(which.max(t)))
    r <- bottom$pointID[bottom$clusterID != i]
    tree <- tree[-r]
  }

  toc(ti)

  return(tree)
}


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
  t0 <- tic() ; gc()
  #step = 0.2; cl_dist = 0.1; max_d = 0.1; apex = 0.005

  j    <- which.min(tree$Z)
  tx   <- tree$X[j]
  ty   <- tree$Y[j]
  tz   <- tree$Z[j]

  tree <- shift(tree, tx, ty, tz)   # Move to origin for numerical stability
  tree <- filter_tree(tree)
  tree <- clean_tree_butt(tree)

  qsm  <- qsm_skeleton(tree, step, cl_dist, max_d)
  qsm  <- qsm_architecture(qsm)
  qsm  <- qsm_detect_weird_butt(qsm)

  d    <- estimate_prolongation(tree, qsm)

  qsm  <- qsm_prolongation(qsm, d)
  qsm  <- qsm_radius(qsm, tree, R0, tip_radius = apex)
  qsm  <- qsm_volume(qsm)
  qsm  <- shift(qsm, -tx, -ty, -tz)

  toc(t0, space = "")

  return(qsm)
}

clean_tree_butt = function(tree)
{
  cat("Cleaning tree butt...") ; ti = tic()

  tree@data$pointID <- 1:lidR::npoints(tree)
  bottom <- tree[tree$Z < min(tree$Z) + 0.25]
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

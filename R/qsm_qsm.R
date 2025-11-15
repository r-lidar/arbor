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

  qsm <- qsm_build_skeleton(tree, step, cl_dist, max_d, verbose)
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

  R0  <- find_root_radius(tree, qsm)
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
  cat("Cleaning tree butt\n") ; ti = tic()

  tree@data$pointID <- 1:lidR::npoints(tree)
  bottom <- tree[tree$Z < 0.25]
  bottom <- lidR::connected_components(bottom, 0.05, 10, connectivity = 26)

  if (length(unique(bottom$clusterID)) > 1)
  {
    cat("\033[33m  Multiple clusters at the bottom of the tree detected. Automatic cleaning triggered.\033[0m\n")

    t <- table(bottom$clusterID)
    i <- as.numeric(names(which.max(t)))
    r <- bottom$pointID[bottom$clusterID != i]
    tree <- tree[-r]
  }

  toc(ti)

  return(tree)
}

qsm_remove_disconnected_branches <- function(dt)
{
  # Ensure data.table
  if (!data.table::is.data.table(dt)) dt <- data.table::as.data.table(dt)

  # Validate required columns
  required_cols <- c("cyl_ID", "parent_ID", "axis_ID")
  if (!all(required_cols %in% names(dt))) {
    stop("Input must contain columns: cyl_ID, parent_ID, axis_ID")
  }

  # Step 1: start with cylinders on the main axis
  keep <- dt[axis_ID == 1, cyl_ID]
  new <- keep

  # Step 2: recursively find all descendants
  repeat {
    children <- dt[parent_ID %in% new, cyl_ID]
    new <- setdiff(children, keep)
    if (length(new) == 0) break
    keep <- c(keep, new)
  }

  # Step 3: keep only connected cylinders
  dt[cyl_ID %in% keep]
}

qsm_detect_weird_butt = function(qsm)
{
  cat("Validating butt architecture \n")  ; ti = tic()

  qsm$angle <- with(qsm,{
    dx <- endX - startX
    dy <- endY - startY
    dz <- endZ - startZ
    acos(dz / sqrt(dx^2 + dy^2 + dz^2)) * 180 / pi
  })

  main   <- qsm[axis_ID == 1]
  angles <- main$angle
  thresh <- 50
  window <- 4  # number of consecutive values required below threshold
  i      <- 1
  while (i <= length(angles))
  {
    if (all(angles[i:min(i+window-1, length(angles))] < thresh)) break
    i <- i + 1
  }

  if (i > 1)
  {
    cat("\033[33m  Detection of weird tree butt. Automatic fix triggered.\033[0m\n")

    main <- qsm[axis_ID == 1]
    rm   <- main[1:i]
    qsm  <- qsm[!cyl_ID %in% rm$cyl_ID]
    qsm  <- qsm_remove_disconnected_branches(qsm)
    j <- which(qsm$axis_ID == 1)[1]
    qsm[j, parent_ID := 0]
  }

  toc(ti)

  return(qsm)
}


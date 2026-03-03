clean_tree_butt = function(tree)
{
  logger("Cleaning tree butt")
  data = qsm_clean_tree_butt_cpp(tree@data);
  data.table::setDT(data)
  tree@data = data
  tree = lidR::las_update(tree);
  return(tree)

  # Old R code
  # tree@data$pointID <- 1:lidR::npoints(tree)
  # bottom <- tree[tree$Z < min(tree$Z) + 0.25]
  # bottom <- connected_components(bottom, 0.05, 10, connectivity = 26)
  #
  # if (length(unique(bottom$clusterID)) > 1)
  # {
  #   logger("Multiple clusters at the bottom of the tree detected. Automatic cleaning triggered.", level = "WARN")
  #
  #   t <- table(bottom$clusterID)
  #   i <- as.numeric(names(which.max(t)))
  #   r <- bottom$pointID[bottom$clusterID != i]
  #   tree <- tree[-r]
  # }
  #
  # return(tree)
}

qsm_detect_weird_butt = function(qsm)
{
  axis_ID <- cyl_ID <- parent_ID <- NULL

  logger("Validating butt architecture")

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
    logger("Detection of weird tree butt. Automatic fix triggered.", level = "WARN")

    main <- qsm[axis_ID == 1]
    rm   <- main[1:i]
    qsm  <- qsm[!cyl_ID %in% rm$cyl_ID]
    qsm  <- qsm_remove_disconnected_branches(qsm)
    j <- which(qsm$axis_ID == 1)[1]
    qsm[j, parent_ID := 0]
  }

  return(qsm)
}

qsm_remove_disconnected_branches <- function(dt)
{
  axis_ID <- cyl_ID <- parent_ID <- NULL

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

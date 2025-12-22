make_master_seed = function(las)
{
  x <- mean(las$X)
  y <- mean(las$Y)
  z <- min(las$Z)-1
  seed <- data.frame(X= x, Y = y, Z = z, anisotropy = 1, pointID = 0)
  lidR::quantize(seed[["X"]], 0.01, las@header[["X offset"]])
  lidR::quantize(seed[["Y"]], 0.01, las@header[["Y offset"]])
  lidR::quantize(seed[["Z"]], 0.01, las@header[["Z offset"]])
  header <- rlas::header_create(seed)
  seed   <- suppressWarnings(lidR::LAS(seed, header))
  seed
}

evaluate_penalty = function(params)
{
  penalty = params$path_finder$angle_penalty(0:180)
  if (length(penalty) != 181) stop("Invalid penalty function")
  params$path_finder$penalty = penalty
  params
}

connected_components = function(las, res, min_pts, name = "clusterID", connectivity = 6)
{
  .N <- N <- clusterID <- NULL

  allowed <- c(6L, 18L, 26L)
  if (!connectivity %in% allowed) stop(sprintf("Invalid connectivity: %d. Allowed values are %s", connectivity, paste(allowed, collapse = ", ")))

  u = C_connected_component(las@data, res, connectivity)
  las = add_lasattribute(las, u, name, "connected component ID")
  grp = las@data[, .N, by = clusterID]
  grp = grp[N < min_pts]
  invalid = las@data[[name]] %in% grp$clusterID
  las@data[[name]][invalid] = 0L
  return(las)
}

sor = function(las, k, m)
{
  noise = C_sor(las@data, k, m, lidR::get_lidr_threads())

  if ("Classification" %in% names(las))
  {
    new_classes <- las@data[["Classification"]]
    new_classes[new_classes == lidR::LASNOISE] <- lidR::LASUNCLASSIFIED
  }
  else
  {
    new_classes <- rep(lidR::LASUNCLASSIFIED, lidR::npoints(las))
  }

  new_classes[noise] <- lidR::LASNOISE
  las@data[["Classification"]] <- new_classes
  return(las)
}


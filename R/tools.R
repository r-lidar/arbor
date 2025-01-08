seed_from_dtm = function(dtm)
{
  seeds = terra::rast(terra::ext(dtm), res = 0.25)
  seeds = terra::resample(dtm, seeds)
  seeds = as.data.frame(seeds, xy = T)
  seeds = data.table::as.data.table(seeds)
  names(seeds) = c("X", "Y", "Z")
  seeds
}

compute_network = function(data, query, k = 5)
{
  if (missing(query))
  {
    nn = lidR::knn(data, k = k)
    n = npoints(data)
  }
  else
  {
    nn = lidR::knnx(data, query, k = k)
    n = npoints(query)
  }

  from <- rep(1:n, each = k)
  to <- as.vector(t(nn$nn.index))
  cost <- as.vector(t(nn$nn.dist))
  edges <- data.frame(from, to, cost)
  edges
}

expand_treeid_to_neighbors = function(unclustered, clustered, max_gap = 0.5, z_factor = 1)
{
  ID = "treeID"
  clustered$Z = clustered$Z*z_factor
  unclustered$Z = unclustered$Z*z_factor
  seed_to_dense_lookup <- knnx(clustered, unclustered, k = 1)
  full_tree_id_vector <- clustered@data[[ID]][seed_to_dense_lookup[[1]]]
  distances = seed_to_dense_lookup[[2]]
  full_tree_id_vector[distances > max_gap] <- NA_integer_
  unclustered <- lidR::add_lasattribute(unclustered, full_tree_id_vector, name = "treeID", desc = "tree ID")
  unclustered$Z = unclustered$Z/z_factor
  return(unclustered)
}

tic = function()
{
  Sys.time()
}

toc = function(t0, space = "  ")
{
  units = "secs"
  tf = Sys.time()
  dt = difftime(tf, t0, units = units)

  # Check if time is greater than 120 seconds and switch to minutes
  if (dt > 120) {
    dt = dt / 60  # Convert to minutes
    units = "mins"
  }

  cat(paste0(space, "Done in ", round(dt, 1), " ", units,  "\n"))
}

free = function(...)
{
  object_names <- as.character(substitute(list(...)))[-1L]

  bytes = 0
  for (obj in list(...)) {
    bytes = bytes + object.size(obj)
  }
  mb = format(bytes, "MB")

  # Remove the specified objects
  for (obj in object_names) {
    if (exists(obj, envir = parent.frame())) {
      rm(list = obj, envir = parent.frame())
    } else {
      warning(paste("Object", obj, "does not exist."))
    }
  }

  gc()

  cat("  Memory freed:", mb, "\n")
}

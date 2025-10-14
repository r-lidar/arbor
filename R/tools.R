seed_from_dtm = function(dtm, res)
{
  seeds = terra::rast(terra::ext(dtm), res = res)
  seeds = terra::resample(dtm, seeds)
  seeds = as.data.frame(seeds, xy = T)
  seeds = data.table::as.data.table(seeds)
  names(seeds) = c("X", "Y", "Z")
  seeds
}


tic = function()
{
  Sys.time()
}

toc = function(t0, space = "  ")
{
  gc()
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
    bytes = bytes + utils::object.size(obj)
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
  #cat("  Memory freed:", mb, "\n")
}

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

is.tls = function(las)
{
  return(lidR::sensor(las) == lidR:::TLSLAS)
}

stop_if_not_tls = function(las)
{
  if (!is.tls(las))
    stop("This point cloud is not flagged as TLS. It has not been read with lidR::readTLSLAS(). Please use the correct read function.")
}

free = function(...)
{
  object_names <- as.character(substitute(list(...)))[-1L]

  bytes = utils::object.size(NULL)
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

  ans = gc()
  return(invisible())
  #cat("  Memory freed:", mb, "\n")
}

cc_export_labels = function(heights, circles, file)
{
  labels1 = heights[,1:3]
  labels1$label = paste0(round(heights$Z, 1), " m")

  labels2 = circles[,1:3]
  labels2$label = paste0(round(circles$radius*2*100, 1), " cm")
  names(labels2) = names(labels1)

  labels = rbind(labels1, labels2)

  circles = sf::st_as_sf(circles, coords = c("center_x", "center_y"))
  circles = sf::st_buffer(circles, circles$radius)
  circles = sf::st_cast(circles, "POINT")

  path = tools::file_path_sans_ext(file)

  flabel = paste0(path, "_labels.txt")
  frings = paste0(path, "_rings.shp")

  data.table::fwrite(labels, flabel)
  sf::st_write(circles["center_z"], frings, append = FALSE)
}

record_entry <- function(las, timing, file_path) {

  # Hardcoded database path
  .DB_PATH <- "/home/jr/Documents/r-lidar inc/arbor/processing_log.csv"

  # --- Input validation ---
  stopifnot(
    "las must be a LAS object"        = inherits(las, "LAS"),
    "timing must be a difftime object" = inherits(timing, "difftime"),
    "file_path must be a character"    = is.character(file_path) && nchar(file_path) > 0
  )

  # --- Extract metrics from the LAS object ---
  n_points <- nrow(las@data)
  area_m2  <- round(sf::st_area(las), 4)          # lidR::area(), in m²

  # --- Normalise timing to seconds (numeric) ---
  timing_sec <- as.numeric(timing, units = "secs")

  # --- Build the new entry ---
  new_entry <- data.frame(
    file_path  = file_path,
    n_points   = n_points,
    area_m2    = area_m2,
    timing_sec = timing_sec,
    date       = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )

  # --- Read existing database (or create an empty one) ---
  if (file.exists(.DB_PATH)) {
    db <- utils::read.csv(.DB_PATH, stringsAsFactors = FALSE)
  } else {
    message("Database not found - creating a new one at: ", .DB_PATH)
    db <- data.frame(
      file_path  = character(),
      n_points   = integer(),
      area_m2    = numeric(),
      timing_sec = numeric(),
      date       = character(),
      stringsAsFactors = FALSE
    )
  }

  # --- Remove any existing entry for this file_path ---
  n_before <- nrow(db)
  db <- db[db$file_path != file_path, ]
  n_removed <- n_before - nrow(db)
  if (n_removed > 0) {
    message(sprintf("Removed %d existing entry/entries for: %s", n_removed, file_path))
  }

  # --- Append the new entry ---
  db <- rbind(db, new_entry)

  # --- Save back to disk ---
  utils::write.csv(db, file = .DB_PATH, row.names = FALSE)
  message(sprintf(
    "Entry recorded - file: %s | points: %d | area: %.2f m\u00b2 | timing: %.2f s | date: %s",
    file_path, n_points, area_m2, timing_sec, new_entry$date
  ))

  invisible(db)
}



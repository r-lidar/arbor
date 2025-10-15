qsm_topology = function(dt)
{
  dt <- data.table::as.data.table(dt)
  dt$cyl_ID = 1:nrow(dt)

  # Add a column to store parent_ID
  dt[, parent_ID := 0L]

  # Create a key from end coordinates for fast lookup (rounded to avoid floating-point issues)
  coord_key <- function(x, y, z, digits = 6)
  {
    paste0(round(x, digits), "_", round(y, digits), "_", round(z, digits))
  }

  # Create lookup table using segment ends
  dt[, end_key := coord_key(endX, endY, endZ)]
  end_lookup <- dt[, .(cyl_ID, end_key)]
  data.table::setkey(end_lookup, end_key)

  # Match start points to ends of other segments
  dt[, start_key := coord_key(startX, startY, startZ)]
  dt[, parent_ID := end_lookup[.SD, on = "end_key==start_key", cyl_ID]]

  # Replace NAs with 0 (i.e., root segments)
  dt[is.na(parent_ID), parent_ID := 0L]

  # Optional: remove keys
  dt[, c("start_key", "end_key") := NULL]
  dt
}

qsm_topology = function(dt)
{
  dt <- data.table::as.data.table(dt)
  dt[, cyl_ID := .I]

  # Add parent_ID column
  dt[, parent_ID := 0L]

  # Function to create coordinate keys (rounded to avoid floating-point issues)
  coord_key <- function(x, y, z, digits = 6)
  {
    paste0(round(x, digits), "_", round(y, digits), "_", round(z, digits))
  }

  # Create lookup table using segment ends
  dt[, end_key := coord_key(endX, endY, endZ)]
  end_lookup <- dt[, .(cyl_ID, end_key)]
  data.table::setkey(end_lookup, end_key)

  # Create start keys
  dt[, start_key := coord_key(startX, startY, startZ)]

  # Lookup: find the cyl_ID in end_lookup whose end_key matches this start_key
  match_parent <- end_lookup[dt, on = .(end_key = start_key), mult = "first"]

  # Extract the cyl_IDs from the joined table (aligned with dt)
  dt[, parent_ID := match_parent$cyl_ID]

  # Replace NAs with 0 (root segments)
  dt[is.na(parent_ID), parent_ID := 0L]

  # Clean up
  dt[, c("start_key", "end_key") := NULL]

  dt
}

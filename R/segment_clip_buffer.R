#' Clip Trees Using a Buffered Convex Hull
#'
#' This function removes trees located near the edges of a point cloud by clipping them
#' using a buffered convex hull. It first computes the convex hull of the input LAS object
#' and shrinks it by the specified buffer distance. Only trees with seed points inside
#' this buffered region are retained.
#'
#' @param las A LAS object from lidR containing segmented trees.
#' @param buffer Numeric value (in meters). The distance by which the convex hull is shrunk
#'   before filtering trees. Default is -5 (removes trees within 5 meters of the boundary).
#' @return A filtered LAS object containing only trees whose seeds fall within the buffered region.
#' @export
clip_buffer = function(las, buffer = -5)
{
  # Avoid NOTES
  X <- Y <- Z <- treeID <- NULL

  # Make root robust to small groups
  root <- function(x, y, z)
  {
    i <- order(z)
    n <- min(100L, length(i))            # Use available points (<=100)
    if (n < 10L) {                       # Guard against tiny groups
      return(list(X = NA_real_, Y = NA_real_, Z = NA_real_))
    }
    x <- x[i][seq_len(n)]
    y <- y[i][seq_len(n)]
    z <- z[i][seq_len(n)]
    x <- mean(x)
    y <- mean(y)
    z <- mean(z)
    return(list(X = x, Y = y, Z = z))
  }

  # Compute roots per tree (drop unsegmented)
  roots <- las@data[!is.na(treeID), root(X, Y, Z), by = treeID]
  roots <- roots[is.finite(X) & is.finite(Y)]     # Drop invalid coords
  if (!nrow(roots)) return(las[0])                # Early exit if nothing valid

  # Build buffered convex hull
  bb <- suppressWarnings(sf::st_convex_hull(las))
  bb <- suppressWarnings(sf::st_buffer(bb, dist = buffer))
  if (any(sf::st_is_empty(bb))) return(las[0])    # Early exit if buffer kills hull

  # Convert roots to sf (safe: no NA)
  seeds <- sf::st_as_sf(roots, coords = c("X", "Y"), crs = sf::st_crs(las))
  valid_seeds <- suppressWarnings(sf::st_filter(seeds, bb))
  if (!nrow(valid_seeds)) return(las[0])

  # Keep trees whose seeds are inside
  las[las$treeID %in% valid_seeds$treeID]
}

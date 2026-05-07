#' Remove Small Trees and Clean Understory
#'
#' This function removes small trees that are likely to be poorly segmented and filters out
#' low understory vegetation based on height. The function does not removes the points. Instead the
#' tree IDs become negative. Negative tree IDs are never processed by \link{qsf} can easily be reverted.
#' See the [Arbor book](<placeholder>) for mode details.
#'
#' @param las A LAS object from lidR.
#' @param max_height Trees with a height less than this threshold (in meters) will be removed.
#' @export
remove_small_trees = function(las, max_height = 2)
{
  treeID <- hag <- hag_max <- hag_min <- NULL

  attributes <- names(las)
  stopifnot("treeID" %in% attributes)
  stopifnot("hag" %in% attributes)

  ans   <- las@data[!is.na(treeID), list(hag_max = max(hag), hag_min = min(hag)), by = treeID]
  ans   <- ans[hag_max > max_height & hag_min < max_height]

  # Keep trees whose seeds are inside
  las@data$treeID = data.table::copy(las@data$treeID)
  las@data[treeID < 0, treeID := treeID * -1] # revert previous
  las@data[!treeID %in% ans$treeID, treeID := treeID * -1]
  las
}

keep_small_trees = function(las, max_height = 2)
{
  treeID <- hag <- hag_max <- hag_min <- NULL

  attributes <- names(las)
  stopifnot("treeID" %in% attributes)
  stopifnot("hag" %in% attributes)

  ans   <- las@data[!is.na(treeID), list(hag_max = max(hag), hag_min = min(hag)), by = treeID]
  ans   <- ans[hag_max < max_height]

  trees <- lidR::filter_poi(las, treeID %in% ans$treeID & !is.na(treeID))

  trees
}


#' Clip Trees Using a Buffer
#'
#' Removes trees located near the edges of a point cloud by clipping the tree that are beyond the
#' limit of the polygon. It first computes the convex hull of the input LAS object
#' and shrinks it by the specified buffer distance. Only trees with seed points inside
#' this buffered region are retained. The function does not removes the points. Instead the tree IDs
#' of trees whose seeds fall outside the region of interest become negative. Negative tree IDs are
#' never processed by \link{qsf} and buffering can easily be reverted. See the [Arbor book](<placeholder>)
#' for mode details.
#'
#' @param las A LAS object from lidR containing segmented trees.
#' @param seeds A LAS object. The seeds from \link{find_seeds}. If missing an internal routine
#' will estimate the position of the trees based on their lowest points.
#' @param buffer Numeric value (in meters). The distance by which the convex hull is shrunk
#'   before filtering trees. Default is -5 (removes trees within 5 meters of the boundary).
#'   Can also be a sf POLYGON object to clip a more complex polygon.
#' @export
#' @md
#' @importFrom data.table :=
remove_buffer = function(las, seeds, buffer = -5)
{
  if (!"treeID" %in% names(las))   stop("Input point cloud must have an attribute 'treeID'")
  if (!"treeID" %in% names(seeds)) stop("Input seeds must have an attribute 'treeID'")

  # Avoid NOTES
  X <- Y <- Z <- treeID <- NULL

  # Compute roots per tree (drop unsegmented)
  if (missing(seeds))
  {
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

    roots <- las@data[!is.na(treeID), root(X, Y, Z), by = treeID]
    roots <- roots[is.finite(X) & is.finite(Y)]     # Drop invalid coords
    if (!nrow(roots)) return(las[0])                # Early exit if nothing valid
  }
  else
  {
    # Make root robust to small groups
    root <- function(x, y, z)
    {
      x <- mean(x)
      y <- mean(y)
      z <- mean(z)
      return(list(X = x, Y = y, Z = z))
    }

    roots <- seeds@data[!is.na(treeID), root(X, Y, Z), by = treeID]
    roots <- roots[is.finite(X) & is.finite(Y)]     # Drop invalid coords
    if (!nrow(roots)) return(las[0])                # Early exit if nothing valid
  }

  # Build buffered convex hull
  if (is.numeric(buffer))
  {
    bb <- suppressWarnings(sf::st_convex_hull(las))
    bb <- suppressWarnings(sf::st_buffer(bb, dist = buffer))
    if (any(sf::st_is_empty(bb))) return(las[0])    # Early exit if buffer kills hull
  }

  if (inherits(buffer, "sf") || inherits(buffer, "sfc"))
  {
    if (!any(sf::st_geometry_type(buffer) %in% c("POLYGON", "MULTIPOLYGON")))
      stop("buffer must be a polygon or multipolygon")

    bb <- buffer
  }

  # Convert roots to sf (safe: no NA)
  seeds <- sf::st_as_sf(roots, coords = c("X", "Y"), crs = sf::st_crs(las))
  valid_seeds <- suppressWarnings(sf::st_filter(seeds, bb))
  if (!nrow(valid_seeds)) return(las[0])

  # Keep trees whose seeds are inside
  las@data$treeID = data.table::copy(las@data$treeID)
  las@data[treeID < 0, treeID := treeID * -1] # revert previous
  las@data[!treeID %in% valid_seeds$treeID, treeID := treeID * -1]
  return(las)
}


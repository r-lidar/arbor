#' Clip Trees Using a Buffered Convex Hull
#'
#' This function removes trees located near the edges of a point cloud by clipping them
#' using a buffered convex hull. It first computes the convex hull of the input LAS object
#' and shrinks it by the specified buffer distance. Only trees with seed points inside
#' this buffered region are retained.
#'
#' @param las A LAS object from lidR containing segmented trees.
#' @param circles
#' @param buffer Numeric value (in meters). The distance by which the convex hull is shrunk
#'   before filtering trees. Default is -5 (removes trees within 5 meters of the boundary).
#' @return A filtered LAS object containing only trees whose seeds fall within the buffered region.
#' @export
clip_buffer = function(las, circles, buffer = -5)
{
  bb = sf::st_convex_hull(las)
  bb = sf::st_buffer(bb, dist = buffer)

  seeds = sf::st_as_sf(circles, coords = c("center_x", "center_y"))
  valid_seeds = sf::st_filter(seeds, bb)

  valid_trees = las[las$treeID %in% valid_seeds$treeID]
  valid_trees
}

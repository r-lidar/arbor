#' @export
clip_buffer = function(las, seeds, buffer = -5)
{
  bb = sf::st_convex_hull(las)
  bb = sf::st_buffer(bb, dist = buffer)
  valid_seeds = sf::st_filter(seeds, bb)

  valid_trees = las[las$treeID %in% valid_seeds$treeID]
  valid_trees
}

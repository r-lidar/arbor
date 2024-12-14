#' @export
clip_buffer = function(las, seeds, buffer = -5)
{
  bb = sf::st_bbox(las)
  bb = terra::ext(bb)
  bb = bb-5
  bb = sf::st_bbox(bb)
  valid_seeds = sf::st_crop(seeds, bb)

  valid_trees = las[las$treeID %in% valid_seeds$treeID]
  valid_trees
}

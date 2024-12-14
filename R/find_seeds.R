#' @export
find_seeds = function(las, slice_seeds_at = c(0.5, 0.8))
{
  seed = lidR::filter_poi(las, hag > slice_seeds_at[1], hag < slice_seeds_at[2])
  seed = lidR::filter_poi(seed, foliage == FALSE)

  #plot(seed)

  seed = lidR::connected_components(seed, 0.05, 10)
  seed = lidR::filter_poi(seed, clusterID != 0)

  seeds = seed@data[, .(X = round(mean(X),3), Y = round(mean(Y),3), Z = round(min(Z), 3)), by = clusterID]
  sfseeds = sf::st_as_sf(seeds, coords = c("X", "Y", "Z"))
  names(sfseeds)[1] = "treeID"

  return(sfseeds)
}

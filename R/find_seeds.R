#' @export
find_seeds = function(las, slice_seeds_at = c(0.5, 0.8))
{
  seed = lidR::filter_poi(las, hag > slice_seeds_at[1], hag < slice_seeds_at[2], foliage == FALSE)

  seed = lidR::connected_components(seed, 0.05, 10)
  seed = lidR::filter_poi(seed, clusterID != 0)

  f = function(x,y,z)
  {
    Z = c(round(min(z), 3), round(max(z), 3))
    dZ = diff(Z)

    bottom = z < Z[2] - dZ/2
    top = z > Z[1] + dZ/2

    X = c(round(median(x[bottom]), 3), round(median(x[top]), 3))
    Y = c(round(median(y[bottom]), 3), round(median(y[top]), 3))

    return(list(X = X,Y = Y,Z = Z))
  }

  seeds = seed@data[, f(X,Y,Z), by = clusterID]

  sfseeds = sf::st_as_sf(seeds, coords = c("X", "Y", "Z"))
  names(sfseeds)[1] = "treeID"

  return(sfseeds)
}

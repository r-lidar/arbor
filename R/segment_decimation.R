#' Homogenization of the Point Cloud
#'
#' Homogenization of the point cloud using a hybrid approach that includes
#' Barycentric Voxel Decimation and reinjection of random points. See the [Arbor book].
#'
#' @param las LAS object from lidR
#' @param res Voxel resolution.
#' @export
hybrid_homogeneization = function(las, res = 0.02)
{
  stop_if_not_tls(las)
  decimated <- TRUE
  keep <- C_homogeneization(las@data, res)
  las <- las[keep]
  free(keep)
  return(las)
}

barycentric_decimation = function(las, res)
{
  decimated <- TRUE
  keep <- C_homogeneization(las@data, res, hybrid = FALSE)
  las <- las[keep]
  free(keep)
  return(las)
}

#' Barycentric pre-decimation of a point cloud
#'
#' This function performs a barycentric decimation of a LiDAR point cloud by
#' dividing the space into voxels of a specified resolution and keeping only
#' the point closest to the barycenter of each voxel. The resulting LAS object
#' includes a logical column indicating which points are retained.
#'
#' @param las A LAS object from lidR.
#' @param params list See \link{parameters}.
#' @export
barycentric_predecimation = function(las, params = default_arbor_parameters)
{
  res <- params$decimation$barycentric_predecimation_resolution
  keep <- C_voxel_barycenter_decimate(las@data$X, las@data$Y, las@data$Z, res)
  las@data$decimated <- keep
  free(keep)
  return(las)
}

barycentric_decimation = function(las, res)
{
  decimated <- TRUE
  keep <- C_voxel_barycenter_decimate(las@data$X, las@data$Y, las@data$Z, res)
  las@data$decimated <- keep
  las <- lidR::filter_poi(las, decimated == TRUE)
  las@data$decimated = NULL
  free(keep)
  return(las)
}
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
  decimated <- TRUE
  keep <- C_voxel_barycenter_decimate(las@data$X, las@data$Y, las@data$Z, res)
  n = sum(keep)
  rm = which(!keep)
  i = sample(rm, n*0.1)
  keep[i] = TRUE
  las = las[keep]
  free(keep)
  return(las)
}


get_barycentric_predecimation <- function(las, params = default_arbor_parameters)
{
  decimated <- NULL

  if (!"decimated" %in% names(las))
    las <- barycentric_predecimation(las, params)

  dec <- lidR::filter_poi(las, decimated == TRUE)
  dec
}

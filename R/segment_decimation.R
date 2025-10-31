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
barycentric_predecimation = function(las, params = default_parameters)
{
  res <- params$decimation$barycentric_predecimation_resolution
  id <- lidR:::C_voxel_id(las, res)
  keep <- C_voxel_barycenter_decimate(las@data$X, las@data$Y, las@data$Z, id)
  las@data$decimated <- keep
  free(id, keep)
  return(las)
}

barycentric_decimation = function(las, res)
{
  decimated <- TRUE
  id <- lidR:::C_voxel_id(las, res)
  keep <- C_voxel_barycenter_decimate(las@data$X, las@data$Y, las@data$Z, id)
  las@data$decimated <- keep
  free(id, keep)
  las <- lidR::filter_poi(las, decimated == TRUE)
  return(las)
}

get_barycentric_predecimation <- function(las, params = default_parameters)
{
  decimated <- NULL

  if (!"decimated" %in% names(las))
    dec <- barycentric_predecimation(las, params)

  dec <- lidR::filter_poi(las, decimated == TRUE)
  dec
}

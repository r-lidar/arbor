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
  logger(sprintf("Barycentric predecimation"))
  res <- params$decimation$barycentric_predecimation_resolution
  keep <- C_homogeneization(las@data, res, hybrid = FALSE)
  las@data$decimated <- keep
  logger(sprintf("Retention %d points (%.1f%%)", sum(keep), sum(keep) / length(keep) * 100))
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
  logger(sprintf("Hybrid homogeneization"))
  decimated <- TRUE
  keep <- C_homogeneization(las@data, res)
  las = las[keep]
  logger(sprintf("Retention %d points (%.1f%%)", sum(keep), sum(keep) / length(keep) * 100))
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

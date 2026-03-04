#' Classify ground points
#'
#' This function classifies ground points in a \code{LAS} object using a custom pipeline optimized for MLS/TLS
#' The detected ground points are thinned, denoised to produce a cleaner and more robust ground classification.
#'
#' @param las A \code{LAS} object.
#'
#' @return A \code{LAS} object with updated ground classification.
#'
#' @export
arbor_ground = function(las)
{
  las@data$pointID = 1:lidR::npoints(las)
  algo   <- lidR::csf(rigidness = 1, class_threshold = 0.05, cloth_resolution = 0.1)
  las    <- lidR::classify_ground(las, algo, last_returns = FALSE)
  ground <- lidR::filter_poi(las, Classification == LASGROUND)
  ground <- lidR::decimate_points(ground, lidR::lowest(0.25))
  ground <- lidR::classify_noise(ground, lidR::sor(k = 10, m = 2))
  ground <- lidR::remove_noise(ground)
  ground$Classification <- lidR::LASGROUND
  las@data$Classification[las@data$Classification == lidR::LASGROUND] = lidR::LASUNCLASSIFIED
  las@data$Classification[ground@data$pointID] = lidR::LASGROUND
  las@data$pointID = NULL
  return(las)
}

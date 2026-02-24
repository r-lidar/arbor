#' Compute Local Wood Likelihood
#'
#' This function calculates the local anisotropy of each point in a LAS object using
#' a k-nearest neighbors (KNN) approach. Anisotropy describes the degree of directional
#' variation in the local neighborhood of a point. The number of nearest neighbors used to compute
#' anisotropy depends on point density and sensor accuracy. For datasets with ~15,000–20,000 points/m² and
#' sensor accuracy of 2–4 cm, `k = 75` is suitable. For highly accurate TLS data (millimeter precision)
#' `k` can be reduced. If the density is very high (e.g., 40,000 points/m²), `k` should be increased.
#' The key consideration is that `k` must be large enough to capture local geometric structures but
#' small enough to remain "local." If k = 0 it uses an auto-adaptive value (experimental)
#'
#' It returns a new LAS object with an extra attribute \code{pwood}. This allows any user to use
#' alternative wood likelihood function before to use the function \link{segment_semantic}.
#'
#' @param las A LAS object from lidR.
#' @param params list See \link{parameters}.
#' @export
wood_likelihood = function(las, params = default_arbor_parameters)
{
  logger("Wood likelihood start")
  t0 <- tic()
  k <- params$woodlikelihood$k
  anisotropy <- C_anisotropy(las@data, k, lidR::get_lidr_threads())
  las <- lidR::add_lasattribute_manual(las, anisotropy, "pwood", "wood likelyhood", "float")
  free(anisotropy)
  logger("Wood likelihood completed")
  return(las)
}

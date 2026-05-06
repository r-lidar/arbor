#' Compute Local Wood Likelihood
#'
#' This function calculates the local anisotropy of each point in a LAS object using
#' a k-nearest neighbors (KNN) approach. See the [Arbor book](<placeholder>) for mode details.
#'
#' @param las A LAS object from lidR.
#' @param params list See \link{parameters}.
#' @export
wood_likelihood = function(las, params = arbor_parameters_default)
{
  stop_if_not_tls(las)
  k <- params$woodlikelihood$k
  anisotropy <- C_anisotropy(las@data, k)
  las <- lidR::add_lasattribute_manual(las, anisotropy, "pwood", "wood likelyhood", "float")
  free(anisotropy)
  return(las)
}

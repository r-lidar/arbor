#' Compute Local Anisotropy
#'
#' This function calculates the local anisotropy of each point in a LAS object using
#' a k-nearest neighbors (KNN) approach. Anisotropy describes the degree of directional
#' variation in the local neighborhood of a point.
#'
#' @param las A LAS object from lidR.
#' @param k Integer. The number of nearest neighbors used to compute anisotropy. The optimal value
#' depends on point density and sensor accuracy. For datasets with ~15,000–20,000 points/m² and
#' sensor accuracy of 2–4 cm, `k = 50` is suitable. For highly accurate TLS data (millimeter precision)
#' `k` can be reduced. If the density is very high (e.g., 40,000 points/m²), `k` should be increased.
#' The key consideration is that `k` must be large enough to capture local geometric structures but
#' small enough to remain "local." The default is 0 meaning that is uses an auto-adaptive value.
#' @export
compute_anisotropy = function(las, k = 0)
{
  if (k  <= 0)
  {
    area <- as.numeric(lidR::st_area(las))
    count <- lidR::npoints(las)
    constant = 50/15000
    k = round(((count/area) * constant), 0)
    cat('Auto-adaptive k =', k, "\n")
  }

  t0 = tic()
  eigen = lidR::point_eigenvalues(las, k = k, coeffs = FALSE)
  anisotropy <- (eigen[["eigen_largest"]] - eigen[["eigen_smallest"]]) / eigen[["eigen_largest"]]
  las = lidR::add_lasattribute_manual(las, anisotropy, "anisotropy", "anisotropy", "float")
  free(eigen)
  toc(t0)
  return(las)
}

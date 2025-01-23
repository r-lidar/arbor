#' Compute the local anisotropy
#'
#' Compute the local anisotropy of each point using a knn neighboorhood.
#'
#' @param las a LAS object from lidR
#' @param k integer. k-nearest neighboors. k must be increased with the point density and the accuracy
#' of the sensor. k = 50 fits well for 15-20.000 points/m² with a sensor accuracy of 2-4 cm. When
#' using TLS data accurate within a few millimeter k can be reduced. But if the density is e.g. 40.000
#' pts/m² k must be increased. The ideas is that k must be big enough to capture the local geometry and
#' small enough to be 'local'.
#' @export
compute_anisotropy = function(las, k = 50)
{
  t0 = tic()
  eigen = point_eigenvalues(las, k = k, metrics = T)
  las = add_lasattribute_manual(las, eigen$anisotropy, "anisotropy", "anisotropy", "float")
  free(eigen)
  toc(t0)
  return(las)
}

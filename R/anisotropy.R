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

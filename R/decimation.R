#' @export
barycentric_predecimation = function(las, res)
{
  id = lidR:::C_voxel_id(las, res)
  keep = C_voxel_barycenter_decimate(las@data$X, las@data$Y, las@data$Z, id)
  las@data$decimated = keep
  free(id, keep)
  return(las)
}

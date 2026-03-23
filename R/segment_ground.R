#' Ground Segmentation
#'
#' Classifies ground points and computes height above ground.
#'
#' @param las A `LAS` object (from the `lidR` package) containing the point cloud data
#' @param params list See \link{parameters}.
#'
#' @md
#' @export
#' @seealso \link{find_seeds}, \link{segment_semantic}, \link{segment_instance}
segment_ground = function(las, params = default_arbor_parameters)
{
  params <- evaluate_penalty(params)
  if (!"Classification" %in% names(las)) las@data$Classification = 0L
  las@data$hag = 0
  make_ground_points()
  segment_ground_cpp(las@data, params)
  return(las)
}

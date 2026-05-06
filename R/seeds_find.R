#' Find seeds to perform instance segmentation
#'
#' In order to perform instance segmentation with  \link{segment_instance} we need
#' some seeds with reference treeIDs. Thus function finds the seeds. See
#' the [Arbor book](<placeholder>) for mode details.
#'
#' @param las A LAS object from lidR.
#' @param params list See \link{parameters}.
#' @export
#' @md
find_seeds <- function(las, params)
{
  stop_if_not_tls(las)
  cloud = las@data ; cloud[["treeID"]] = -1L # Allocate memory because C++ assumes memory allocated
  seeds = find_seeds_cpp(cloud, params)
  seeds = suppressWarnings(lidR::LAS(seeds, lidR::header(las)))
  return(seeds)
}

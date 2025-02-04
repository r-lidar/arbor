#' Remove Small Trees and Clean Low Understory
#'
#' This function removes small trees that are likely to be poorly segmented and filters out
#' low understory vegetation based on height.
#'
#' @param las A LAS object from lidR.
#' @param max_height Trees with a height less than this threshold (in meters) will be removed.
#' @export
clean_small_cluster = function(las, max_heigh = 3)
{
  treeID <- hag <- hag_max <- hag_min

  attributes = names(las)
  stopifnot("treeID" %in% attributes)
  stopifnot("hag" %in% attributes)

  ans = las@data[!is.na(treeID), list(hag_max = max(hag), hag_min = min(hag)), by = treeID]
  ans = ans[hag_max > max_heigh & hag_min < max_heigh]

  trees = lidR::filter_poi(las, treeID %in% ans$treeID & !is.na(treeID))
  trees
}

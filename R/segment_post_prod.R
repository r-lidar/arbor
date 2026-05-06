#' Remove Small Trees and Clean Understory
#'
#' This function removes small trees that are likely to be poorly segmented and filters out
#' low understory vegetation based on height.
#'
#' @param las A LAS object from lidR.
#' @param max_height Trees with a height less than this threshold (in meters) will be removed.
#' @export
remove_small_trees = function(las, max_height = 2)
{
  treeID <- hag <- hag_max <- hag_min <- NULL

  attributes <- names(las)
  stopifnot("treeID" %in% attributes)
  stopifnot("hag" %in% attributes)

  ans   <- las@data[!is.na(treeID), list(hag_max = max(hag), hag_min = min(hag)), by = treeID]
  ans   <- ans[hag_max > max_height & hag_min < max_height]

  las$treeID[!las$treeID %in% ans$treeID] = NA_integer_
  las
}

keep_small_trees = function(las, max_height = 2)
{
  treeID <- hag <- hag_max <- hag_min <- NULL

  attributes <- names(las)
  stopifnot("treeID" %in% attributes)
  stopifnot("hag" %in% attributes)

  ans   <- las@data[!is.na(treeID), list(hag_max = max(hag), hag_min = min(hag)), by = treeID]
  ans   <- ans[hag_max < max_height]

  trees <- lidR::filter_poi(las, treeID %in% ans$treeID & !is.na(treeID))

  trees
}

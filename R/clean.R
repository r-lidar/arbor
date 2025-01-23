#' Remove small trees that are likely badly segmented and clean the low understory
#'
#' @param las LAS object from lidR
#' @param max_height remove every tree instance smaller than this height
#' @export
clean_small_cluster = function(las, max_heigh = 3)
{
  ans = las@data[!is.na(treeID), .(hag_max = max(hag), hag_min = min(hag)), by = treeID]
  ans = ans[hag_max > max_heigh & hag_min < max_heigh]

  trees = filter_poi(las, treeID %in% ans$treeID & !is.na(treeID))
  trees
}

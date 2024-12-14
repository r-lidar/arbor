#' @export
clean_small_cluster = function(las, max_heigh = 3)
{
  ans = las@data[!is.na(treeID), .(hag_max = max(hag), hag_min = min(hag)), by = treeID]
  ans = ans[hag_max > max_heigh & hag_min < max_heigh]

  trees = filter_poi(las, treeID %in% ans$treeID & !is.na(treeID))
  trees
}

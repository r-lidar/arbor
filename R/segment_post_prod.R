#' Post production tool
#'
#' Apply this function to fix some issue in the instance segmentation. If you have no issue but apply
#' this function this should not have any effect. This function look at the wood part of each trees
#' an test with a connected component algorithm if there are more than one big cluster. This can
#' happen when a seed is missing with \link{find_seeds}. In this case a few wood point may be assigned
#' to the wrong tree because this was the best seed to reach. Those points are reasigned to foliage class
#' in order to do not mess-up QSM software
#'
#' @param las LAS object from lidR with semantic and instance segmentation
#' @export
fix_small_isolated_low_clusters = function(las)
{
  treeID <- clusterID <- foliage <- hag <- foliage <- NULL

  las@data$pointID = 1:lidR::npoints(las)
  sub = lidR::filter_poi(las, foliage == FALSE, hag < 3)

  for (id in unique(las$treeID))
  {
    cat("Tree", id)
    tt = lidR::filter_poi(sub, treeID == id)
    if (lidR::is.empty(tt))
    {
      cat("\n")
      next
    }
    tt$Z = tt$Z * 0.1
    tt = lidR::connected_components(tt, 0.05, 200)
    tt$Z = tt$Z * 10

    ids = 1
    n = length(unique(tt$clusterID))
    if (n > 1)
    {
      ids = table(tt$clusterID)
      ids = as.numeric(names(ids[which.max(ids)]))
      cat(" :", n, "clusters found. Smaller cluster(s) re-assigned as foliage")
      pid = tt$pointID[tt$clusterID != ids]
      las$foliage[pid] = TRUE
      tt = lidR::filter_poi(tt, clusterID == ids)
      #plot(tt, color = "clusterID")
      #cat("  ", id, "\n")
      #plot(tt, color = "foliage", pal = c("chocolate4", "darkgreen"))
    }
    else
    {
      cat(" ok")
    }

    cat("\n")
  }

  las@data$pointID = NULL

  return(las)
}

#' Remove Small Trees and Clean Understory
#'
#' This function removes small trees that are likely to be poorly segmented and filters out
#' low understory vegetation based on height.
#'
#' @param las A LAS object from lidR.
#' @param max_height Trees with a height less than this threshold (in meters) will be removed.
#' @export
remove_small_trees = function(las, max_height = 3)
{
  treeID <- hag <- hag_max <- hag_min <- NULL

  attributes <- names(las)
  stopifnot("treeID" %in% attributes)
  stopifnot("hag" %in% attributes)

  ans   <- las@data[!is.na(treeID), list(hag_max = max(hag), hag_min = min(hag)), by = treeID]
  ans   <- ans[hag_max > max_height & hag_min < max_height]

  trees <- lidR::filter_poi(las, treeID %in% ans$treeID & !is.na(treeID))

  trees
}

clean_small_cluster = remove_small_trees


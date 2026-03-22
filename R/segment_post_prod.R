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
fix_small_isolated_low_clusters <- function(las)
{
  if (!interactive())
  {
    options(progressr.enable = TRUE)
  }

  treeID <- clusterID <- foliage <- hag <- NULL

  las@data$pointID <- 1:lidR::npoints(las)
  sub <- lidR::filter_poi(las, foliage == FALSE, hag < 3)

  tree_ids <- unique(las$treeID)
  n_trees <- length(tree_ids)

  progressr::handlers(
    progressr::handler_progress(
      format = ":spin :current/:total [:bar] :percent in :elapsed ETA: :eta",
      width = 50,
      complete = "="
    ))

    progressr::with_progress(
    {
      p <- progressr::progressor(steps = n_trees)

      for (id in tree_ids)
      {
        tt <- lidR::filter_poi(sub, treeID == id)
        if (!lidR::is.empty(tt)) {
          tt$Z <- tt$Z * 0.1
          tt <- connected_components(tt, 0.05, 200)
          tt$Z <- tt$Z * 10

          n <- length(unique(tt$clusterID))
          if (n > 1) {
            ids <- table(tt$clusterID)
            ids <- as.numeric(names(ids[which.max(ids)]))
            pid <- tt$pointID[tt$clusterID != ids]
            las$foliage[pid] <- TRUE
          }
        }

        p()  # update progress
      }
    })

  las@data$pointID <- NULL
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


clean_small_cluster = remove_small_trees


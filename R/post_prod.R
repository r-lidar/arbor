#' @export
fix_splited_trees = function(las)
{
  fit_circle_to_seed = function(id)
  {
    keep = slice$treeID == id
    keep[is.na(keep)] = FALSE
    cl = slice[keep]
    if (npoints(cl) < 10) return(NULL)
    circle = fit_circle(cl, num_iterations = 400, inlier_threshold = 0.02)
    if (!is.null(circle$radius) && circle$radius < 0.25 && circle$angle_range > 90) return(data.frame(X = circle$center_x, Y = circle$center_y, R = circle$radius, id = id))
    else return(NULL)
  }

  zmin = min(las$hag)
  offsets = seq(0,1,0.25)

  for (k in 1:length(offsets))
  {
    slice = filter_poi(las, hag > zmin+offsets[k], hag < zmin+offsets[k]+0.25, foliage == FALSE)

    #x = plot(slice, color = "treeID")

    circles = lapply(unique(slice$treeID), fit_circle_to_seed)
    circles = do.call(rbind, circles)

    #for (i in 1:nrow(circles)) add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], 0)

    circles = sf::st_as_sf(circles, coords = c("X", "Y"))
    circles = sf::st_buffer(circles, circles$R*1.20)

    intersect = sf::st_intersects(circles)
    intersect = Filter(function(x) length(x) > 1, intersect)
    intersect = unique(intersect)

    if (length(intersect) == 0) next

    for (i in intersect)
    {
      ids = circles[i,]$id
      #plot(filter_poi(las, treeID %in% ids), color = "treeID", axis = T)
      #plot(filter_poi(las, treeID %in% ids), color = "foliage")
      cat("Combining trees", ids, "\n")
      las$treeID[las$treeID %in% ids] = ids[1]
    }
  }

  return(las)
}

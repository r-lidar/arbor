find_root_radius = function(tree, qsm, verbose = FALSE)
{
  xyz = tree@data[,1:3]
  zoffset = min(tree$Z)

  centroid = list()
  centroid$centroidX = (qsm$startX + qsm$endX)/2
  centroid$centroidY = (qsm$startY + qsm$endY)/2
  centroid$centroidZ = (qsm$startZ + qsm$endZ)/2
  centroid = as.data.frame(centroid)
  nearest = FNN::knnx.index(data = centroid, query = xyz, algorithm = "kd_tree", k = 1)

  xyz$cyl_ID = qsm$cyl_ID[nearest]

  qsm = qsm[cyl_ID > 0]
  xyz = xyz[cyl_ID > 0]

  data.table::setkey(main_axis, cyl_ID)
  data.table::setkey(xyz, cyl_ID)

  #rgl::points3d(xyz, col = lidR:::set.colors(xyz$cyl_ID, lidR::pastel.colors(200)))

  main_axis = qsm[axis_ID == 1]
  cylinder_id_queue <- main_axis$cyl_ID[main_axis$parent_ID == 0] # start from the root
  height = 0
  radii = c()
  # Distance to root
  while (length(cylinder_id_queue) > 0 & height < 2) #& length(radii) < 10) for hagenia
  {
    currentID <- cylinder_id_queue[[1]]
    currentIndex <- which(main_axis$cyl_ID == currentID)

    height = main_axis$startZ[currentIndex] - zoffset

    childIDs <- main_axis[parent_ID == currentID]$cyl_ID
    cylinder_id_queue <- c(cylinder_id_queue[-1], as.list(childIDs))

    radius = measure_cylinder_radius(xyz, main_axis, currentID)
    if (!is.na(radius))
    {
      radii = c(radii, radius)
    }
  }

  if (length(radii) == 0)
    stop("Failure: unable to measure a diameter")

  radius = median(radii)

  if (verbose) cat("Root radius =", radius, '\n')

  return(radius)

  segment = main_axis[cyl_ID == currentID]
  rgl::points3d(xyz, col = lidR:::set.colors(xyz$cyl_ID, lidR::pastel.colors(150)))
  start <- c(segment$startX, segment$startY, segment$startZ)
  end <- c(segment$endX, segment$endY, segment$endZ)
  m = rbind(start, end)
  cyl = rgl::cylinder3d(m, radius, closed = -1, sides = 16)
  rgl::shade3d(cyl)
}

measure_cylinder_radius = function(pc, qsm, id)
{
  axis = qsm[cyl_ID == id]
  sub = as.matrix(pc[cyl_ID == id, 1:3])

  start = as.numeric(axis[,1:3])
  end = as.numeric(axis[,4:6])

  #rgl::points3d(sub, col = pc$cylID)
  #rgl::segments3d(rbind(start, end), lwd = 3)

  M = compute_rotation_matrix(start, end)
  sub = sub %*% t(M)

  cl = dbscan::dbscan(sub, eps = 0.05, minPts = 1)
  ans = table(cl$cluster)
  res = which.max(ans)
  ids = as.numeric(names(res))
  sub = sub[which(cl$cluster == ids),]
  if (nrow(sub) > 3)
  {
    circle = ransac_circle(sub, num_iterations = 400, inlier_threshold = 0.02)
    radius = circle$radius

    if (FALSE)
    {
      plot(sub, asp = 1, main = paste("Arc =", circle$covered_arc_degree,  "Inlier =", round(circle$percentage_inlier*100), "R =", round(radius,2)))
      symbols(circle$center_x, circle$center_y, circles = circle$radius, inches = FALSE, add = TRUE, fg = "red")
      symbols(circle$center_x, circle$center_y, circles = circle$radius-0.02, lty = 3, inches = FALSE, add = TRUE, fg = "red")
      symbols(circle$center_x, circle$center_y, circles = circle$radius+0.02, lty = 3, inches = FALSE, add = TRUE, fg = "red")
    }
  }
  else
  {
    circle = list(radius = 1000, covered_arc_degree = 0, percentage_inlier = 0)
  }


  if ((circle$covered_arc_degree > 100 & circle$percentage_inlier*100 > 50) | circle$radius < 0.04)
  {
    return(radius)
  }
  else
  {
    return(NA)
  }

  return(radius)
}


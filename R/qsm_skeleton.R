qsm_skeleton = function(tree, step = .2, cl_dist = 0.1, max_d = 0.3)
{
  data <- qsm_layers(tree, step)
  data <- qsm_clusters(data, cl_dist)
  skel <- qsm_nodes(data, max_d)
  skel <- qsm_topology(skel)
  return(skel)
}

qsm_layers = function(tree, step)
{
  logger("Computing layers")

  if (step == "auto") {
    dz = diff(range(tree$Z))
    D  = dz / 100
    D  = max(0.05, D)
  } else {
    D = step
  }

  #data = cpp_compute_layers(as.matrix(tree@data), D)
  data = qsm_layers_cpp(tree@data, D)
  data.table::setDT(data)
  return(data[])
}

#' @importFrom data.table :=
qsm_clusters = function(data, cl_dist)
{
  iter <- cluster <- radius <- X <- Y <- Z <- NULL

  logger("Clustering layers")

  first = TRUE
  for (i in sort(unique(data$iter)))
  {
    if (first)
    {
      # at the first iteration the clustering distance (cl_dist) is the radius of the first layer
      # only one cluster at iteration 1 (the tree base)
      data[iter == i, cluster := 1]
      # compute the radius
      data[iter == i, radius := sqrt((X-mean(X))^2 + (Y-mean(Y))^2 + (Z-mean(Z))^2  )]
      # the average radius is the first clustering distance
      cl_d = mean(data$radius[data$iter == i])
      first = FALSE

      if (FALSE)
      {
        tmp = data[iter == i]
        plot(tmp[, .(X,Y)], asp = 1, cex = 0.25)
        points(mean(tmp$X), mean(tmp$Y), col = "red")
      }
    }
    else
    {
      # at subsequent iterations the clustering distance is the average radius
      # of clustered objects in the curent layer

      # keep points in the curent layer
      in_iter = which(data$iter==i)

      # if there are more than one point in the layer -> cluster it
      if (length(in_iter) >= 2)
      {
        # clustering
        #temp = LAS(data[in_iter,1:3])
        cl = dbscan::dbscan(data[in_iter,1:3], eps = cl_dist, minPts = 1)
        #cl = lidR::connected_components(temp, res = cl_dist, min_pts = 1, connectivity = 26)
        cl_id = cl$cluster
        #cl_id = cl@data$clusterID

        #cl = fastcluster::hclust(stats::dist(data[in_iter,1:3]), method = "single")
        #cl_id = stats::cutree(cl, h = cl_d)
        data[in_iter, cluster := cl_id]
        data[in_iter, radius := sqrt((X-mean(X))^2 + (Y-mean(Y))^2 + (Z-mean(Z))^2  ), by = cluster]  # compute the radius for each cluster
        cl_d = mean(data$radius[data$iter == i])/2           # the new clustering distance is the average radius of all clusters

        if (FALSE)
        {
          tmp = data[in_iter]
          plot(tmp[, 1:2], asp = 1, col = cl_id, cex = 0.25)
        }

        # if the computed clustering distance is too small -> replace by the user defined minimum distance
        #if(cl_d < cl_dist) cl_d = cl_dist
      }
      else
      {
        # if there is only one point, there is only one cluster
        data[in_iter,cluster := 1]
      }
    }
  }

  if (FALSE)
  {
    col = lidR:::set.colors(data$cluster+data$iter, pastel.colors(50000))
    rgl::points3d(data$X, data$Y, data$Z, col = col)
    col = lidR:::set.colors(data$radius, height.colors(500))
    rgl::points3d(data$X, data$Y, data$Z, col = col)
  }

  return(data[])
}

qsm_nodes = function(data, max_d)
{
  logger("Building nodes")
  skel = cpp_build_skeleton(data, max_d)
  data.table::setDT(skel)

  if (FALSE)
  {
    x = plot(tree, bg = "white")
    plot_qsm(skel, add = x)
  }

  return(skel)
}

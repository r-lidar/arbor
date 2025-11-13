qsm_build_skeleton = function(tree, step = .2, cl_dist = 0.1, max_d = 0.3, verbose = FALSE)
{
  pc = tree@data

  if (step == "auto")
  {
    dz = diff(range(tree$Z))
    D  = dz / 100
    D  = max(0.05, D)
  }
  else
  {
    D = step
  }

  # =--------------------------------=
  # Step 1. iteratively compute layers
  # =--------------------------------=

  cat("  Computing layers\n") ; ti = tic()
  data = cpp_compute_layers(as.matrix(pc), D)
  data.table::setDT(data)

  if (FALSE)
  {
    q = quantile(data$iter, probs = 0.95)
    v = data$iter
    #v[v>q] = q+1
    col = lidR:::set.colors(v, lidR::random.colors(2500))
    rgl::points3d(data, col = col)

    col = lidR:::set.colors(v, lidR::height.colors(2500))
    rgl::points3d(data, col = col)

    v = data$dist
    col = lidR:::set.colors(v, lidR::height.colors(250))
    rgl::points3d(data, col = col)
  }

  toc(ti, space =  "      ")

  # =--------------------------------------------------------=
  # Step 2. clustering non connected components in each layer
  # =--------------------------------------------------------=

  cat("  Clustering layers\n")  ; ti = tic()

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

  toc(ti, space =  "      ")

  # =--------------------------=
  # Step 3. Build the skeleton
  # =--------------------------=

  cat("  Building skeleton\n")  ; ti = tic()
  skel = cpp_build_skeleton(data, max_d)

  if (FALSE)
  {
    x = plot(tree, bg = "white")
    plot_qsm(skel, add = x)
  }

  # segments length
  skel = qsm_length(skel)
  skel = skel[skel$length > 0,]

  # Build the coordinates matrix: start1, end1, start2, end2, ...
  # coords <- matrix(t(cbind( skel[, c("startX", "startY", "startZ")],skel[, c("endX",   "endY",   "endZ")])), ncol = 3, byrow = TRUE)
  # rgl::open3d()
  # rgl::segments3d(coords, color = "blue", lwd = 2)

  toc(ti, space =  "      ")

  # =------------------------=
  # Step 4. compute qsm_topology
  # =------------------------=

  cat("  Computing qsm_topology\n")  ; ti = tic()

  skel = qsm_topology(skel)

  if (FALSE)
  {
    x = plot_qsm(skel, cylinder = FALSE)
    passage = filter_poi(tree, passage > 0)
    plot(passage, add = x, size = 3)
  }

  # No root is a bug
  if (sum(skel$parent_ID == 0) == 0)
    stop("Internal error: no root found")

  # Two roots are rare but possible if forking at root
  if (sum(skel$parent_ID == 0) > 1)
  {
    xyz = skel[which(skel$parent_ID == 0)[1],]
    xyz$endX = xyz$startX
    xyz$endY = xyz$startY
    xyz$endZ = xyz$startZ
    xyz$startZ = xyz$startZ - 0.001
    skel = rbind(xyz, skel)
    skel = qsm_topology(skel)
  }

  qsm = data.table::data.table(
    startX = skel$startX,
    startY = skel$startY,
    startZ = skel$startZ,
    endX = skel$endX,
    endY = skel$endY,
    endZ = skel$endZ,
    cyl_ID = skel$cyl_ID,
    parent_ID = skel$parent_ID
  )

  toc(ti, space =  "      ")

  return(qsm)
}

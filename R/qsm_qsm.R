#' Generate a QSM from a single tree point cloud
#'
#' This function processes a tree point cloud to generate a Quantitative Structure Model (QSM) using
#' the bottom to top approach to build the skeleton. Then if fits a radius to each edge of the skeleton
#' either with a ransac circle fitting or a least square approach depending on if ransac successfully
#' fit a good circle. Last it applies an AdTree-like allometric model to correct irrelevant cylinder
#' and compute a more organic tree.
#'
#' @param tree A `LAS` object containing a single tree point cloud with required attributes. Only
#' the point labelled as wood will be used for QSM.
#' @param apex last radius of the branch tips
#'
#' @export
qsm = function(tree, step = 0.2, cl_dist = 0.1, max_d = 0.1, apex = 0.005, power = 1.1, pure_model = FALSE, ..., verbose = FALSE)
{
  t0 = tic()
  #step = 0.2; cl_dist = 0.1; max_d = 0.1; apex = 0.005; power = 1.1; pure_model = FALSE; verbose = FALSE

  tree <- filter_tree(tree)

  d <- if ("hag" %in% names(tree)) min(tree$hag) else 0

  # Move to origin for numerical stability
  j           <- which.min(tree$Z)
  tx          <- tree$X[j]
  ty          <- tree$Y[j]
  tz          <- tree$Z[j]
  tree@data$X <- tree@data$X-tx
  tree@data$Y <- tree@data$Y-ty
  tree@data$Z <- tree@data$Z-tz

  cat("Cleaning tree butt\n") ; ti = tic()

  tree@data$pointID <- 1:lidR::npoints(tree)
  bottom <- tree[tree$Z < 0.25]
  bottom <- lidR::connected_components(bottom, 0.05, 10, connectivity = 26)

  if (length(unique(bottom$clusterID)) > 1)
  {
    cat("\033[33m  Multiple clusters at the bottom of the tree detected. Automatic cleaning triggered.\033[0m\n")
    warning("Multiple clusters at the bottom of the tree detected. Automatic cleaning triggered.")
    t <-table(bottom$clusterID)
    i <- as.numeric(names(which.max(t)))
    r <- bottom$pointID[bottom$clusterID != i]
    tree <- tree[-r]
  }

  toc(ti) ; cat("Building skeleton\n") ; ti = tic()

  qsm <- build_skeleton(tree, step, cl_dist, max_d, verbose)

  toc(ti) ; cat("Building architecture\n")  ; ti = tic()

  qsm = qsm_architecture(qsm)
  qsm = qsm_smooth(qsm, niter = 1)
  qsm = qsm_architecture(qsm)

  toc(ti) ; cat("Validating butt architecture \n")  ; ti = tic()

  i <- detect_weird_butt(qsm)

  if (i > 1)
  {
    cat("\033[33m  Detection of weird tree butt. Automatic fix triggered.\033[0m\n")
    warning("Detection of weird tree butt. Automatic fix triggered.")

    main <- qsm[axis_ID == 1]
    rm   <- main[1:i]
    qsm  <- qsm[!cyl_ID %in% rm$cyl_ID]
    qsm  <- remove_disconnected_branches(qsm)
    j <- which(qsm$axis_ID == 1)[1]
    qsm[j, parent_ID := 0]

    # recompute the prolongation d
    if ("hag" %in% names(tree))
    {
      minZ <- min(qsm$startZ)
      hag  <- tree$hag[tree$Z < minZ]
      d    <- max(hag)
    }
  }

  toc(ti) ; cat("Measuring diameters\n") ; ti = tic()

  R0  <- find_root_radius(tree, qsm, verbose)
  qsm <- qsm_prolongation(qsm, d)

  if (pure_model)
    qsm <- qsm_radius_model(qsm, tree, R0, tip_radius = apex, power = power)
  else
    qsm <- qsm_radius(qsm, tree, R0, tip_radius = apex, power = power)

  toc(ti)

  qsm$startX  <- qsm$startX+tx
  qsm$startY  <- qsm$startY+ty
  qsm$startZ  <- qsm$startZ+tz
  qsm$endX    <- qsm$endX+tx
  qsm$endY    <- qsm$endY+ty
  qsm$endZ    <- qsm$endZ+tz
  tree@data$X <- tree@data$X+tx
  tree@data$Y <- tree@data$Y+ty
  tree@data$Z <- tree@data$Z+tz

  toc(t0, space = "")

  return(qsm)

  #plot_qsm(qsm, color = "branch_order")

  x = plot_qsm(qsm)
  plot(tree, add = x, size = 2, pal = "brown")
}

build_skeleton = function(tree, step = .2, cl_dist = 0.1, max_d = 0.3, verbose = FALSE)
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

measure_cylinder_radius = function(pc, qsm, id)
{
  axis = qsm[.(id)]
  sub = as.matrix(pc[.(id), 1:3])

  start = as.numeric(axis[,1:3])
  end = as.numeric(axis[,4:6])

  #rgl::points3d(sub)
  #rgl::segments3d(rbind(start, end), lwd = 3)

  M = compute_rotation_matrix(start, end)
  sub = sub %*% t(M)

  if (nrow(sub) > 3)
  {
    circle = lidRtls:::ransac_circle(sub, num_iterations = 50, inlier_threshold = 0.015)
    radius = circle$radius

    #plot(sub, asp = 1)
    #symbols(circle$center_x, circle$center_y, circles = circle$radius, inches = FALSE, add = TRUE)
    #symbols(circle$center_x, circle$center_y, circles = circle$radius-0.015, lty = 3, inches = FALSE, add = TRUE, col = "red")
    #symbols(circle$center_x, circle$center_y, circles = circle$radius+0.015, lty = 3, inches = FALSE, add = TRUE, col = "red")
  }
  else
  {
    circle = list(radius = 1000, covered_arc_degree = 0, percentage_inlier = 0)
  }


  if (circle$covered_arc_degree > 100 && circle$percentage_inlier*100 > 0.5)
  {
    return(radius)
  }
  else
  {
    xc = mean(sub[,1])
    yc = mean(sub[,2])
    return(mean(sqrt((sub[,1] - xc)^2 +     (sub[,2] - yc)^2)))
  }

  return(radius)
}

remove_disconnected_branches <- function(dt) {
  # Ensure data.table
  if (!data.table::is.data.table(dt)) dt <- data.table::as.data.table(dt)

  # Validate required columns
  required_cols <- c("cyl_ID", "parent_ID", "axis_ID")
  if (!all(required_cols %in% names(dt))) {
    stop("Input must contain columns: cyl_ID, parent_ID, axis_ID")
  }

  # Step 1: start with cylinders on the main axis
  keep <- dt[axis_ID == 1, cyl_ID]
  new <- keep

  # Step 2: recursively find all descendants
  repeat {
    children <- dt[parent_ID %in% new, cyl_ID]
    new <- setdiff(children, keep)
    if (length(new) == 0) break
    keep <- c(keep, new)
  }

  # Step 3: keep only connected cylinders
  dt[cyl_ID %in% keep]
}

detect_weird_butt = function(qsm)
{
  qsm$angle <- with(qsm,{
    dx <- endX - startX
    dy <- endY - startY
    dz <- endZ - startZ
    acos(dz / sqrt(dx^2 + dy^2 + dz^2)) * 180 / pi
  })

  main   <- qsm[axis_ID == 1]
  angles <- main$angle
  thresh <- 50
  window <- 4  # number of consecutive values required below threshold
  i      <- 1
  while (i <= length(angles))
  {
    if (all(angles[i:min(i+window-1, length(angles))] < thresh)) break
    i <- i + 1
  }

  return(i)
}


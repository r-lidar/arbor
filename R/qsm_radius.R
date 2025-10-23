qsm_radius_adtree = function(qsm, tree, R0, tip_radius = 0.005, power = 1.1)
{
  root = which(qsm$parent_ID == 0)
  w0 = qsm[["subtree_length"]][root]
  wi = qsm[["subtree_length"]]

  # AdTree allometric model
  s = (wi / w0)^power
  s_min = min(s)
  s_max = max(s)
  r1 = tip_radius + (s - s_min) / (s_max - s_min) * (R0 - tip_radius)

  qsm[["radius"]] = r1
  return(qsm)
}

qsm_radius_hagenia = function(qsm, tree, R0, tip_radius = 0.005, power = 1.1)
{
  plot = FALSE

  root = which(qsm$parent_ID == 0)
  w0 = qsm[["subtree_length"]][root]
  wi = qsm[["subtree_length"]]

  # AdTree allometric model
  s = (wi / w0)^power
  s_min = min(s)
  s_max = max(s)
  r1 = tip_radius + (s - s_min) / (s_max - s_min) * (R0 - tip_radius)

  # Actual measurement for each cylindre using RANSAC if it work or falling back
  # to barycenter + average distance
  centroid = list()
  centroid$centroidX = (qsm$startX + qsm$endX)/2
  centroid$centroidY = (qsm$startY + qsm$endY)/2
  centroid$centroidZ = (qsm$startZ + qsm$endZ)/2
  centroid = as.data.frame(centroid)
  nearest = FNN::knnx.index(data = centroid, query = tree@data[,1:3], algorithm = "kd_tree", k = 1)

  xyz = tree@data[, 1:3]
  xyz$cyl_ID = qsm$cyl_ID[nearest]

  if (plot)
  {
    plot_qsm(qsm, cylinder = FALSE)
    rgl::points3d(centroid, size = 3)
    rgl::points3d(xyz, col = lidR:::set.colors(xyz$cyl_ID, pal = lidR::random.colors(nrow(centroid))))
  }

  data.table::setkey(qsm, cyl_ID)
  data.table::setkey(xyz, cyl_ID)

  r2 = numeric(nrow(qsm))
  ransac = logical(nrow(qsm))

  fit_circle_least_squares <- function(xy)
  {
    x <- xy[,1]
    y <- xy[,2]
    A <- cbind(2*x, 2*y, rep(1, length(x)))
    b <- x^2 + y^2
    sol <- qr.solve(A, b)
    xc <- sol[1]
    yc <- sol[2]
    r  <- sqrt(sol[3] + xc^2 + yc^2)
    list(xc = xc, yc = yc, radius = r)
  }

  for (i in 1:nrow(qsm))
  {
    id = qsm$cyl_ID[i]
    axis = qsm[.(id)]
    sub = as.matrix(xyz[.(id), 1:3])

    start = as.numeric(axis[,1:3])
    end = as.numeric(axis[,4:6])

    #rgl::points3d(sub)
    #rgl::segments3d(rbind(start, end), lwd = 3)

    M = compute_rotation_matrix(start, end)
    sub = sub %*% t(M)

    if (nrow(sub) < 5) next

    r = lidRtls:::ransac_circle(sub, num_iterations = 100)

    if (r$covered_arc_degree < 100 | r$percentage_inlier*100 < 50 | r$percentage_inside > 25)
    {
      r2[i] = fit_circle_least_squares(sub[,1:2])$radius
    }
    else
    {
      r2[i] = r$radius
      ransac[i] = TRUE
    }

    #rgl::points3d(sub)
    #lidR::add_circle3d(c(0,0), xc, yc, r, mean(sub[,3]))
  }

  qsm[["radius"]] = r2
  #qsm$ransac = ransac

  if (plot)
  {
    plot_qsm(qsm, cylinder = TRUE)
    rgl::points3d(centroid, size = 3)
    rgl::points3d(xyz, col = lidR:::set.colors(xyz$cyl_ID, pal = lidR::random.colors(nrow(centroid))))
  }

  # Pass through each axis and find the location where measurement strat too be too small to be relevant
  # replace subsequent radii by the allometric model
  for (id in unique(qsm$axis_ID))
  {
    axis = qsm[axis_ID == id]
    if (id == 1) next
    i = which(axis$radius > 0.0 & axis$radius < 0.02)[1]
    if (is.na(i)) next
    w0 = axis$subtree_length[i]
    wi = axis$subtree_length
    s = (wi / w0)^power
    r = s * 0.02 + tip_radius
    r[1:i] = axis$radius[1:i]
    axis$radius = r
    qsm$radius[qsm$axis_ID == id] = r
  }

  # Replace irrelevant radii by the allometric model
  qsm$radius = ifelse(qsm$radius > 0.01 & qsm$radius < r1, qsm$radius, r1)
  qsm$radius = ifelse(qsm$radius > 1.5*r1, r1, qsm$radius)
  qsm$radius = ifelse(is.na(qsm$radius), r1, qsm$radius)

  # Smooth everything
  smooth_radius <- function(x, window = 5, fun = median)
  {
    if (!is.numeric(x) || !is.vector(x)) stop("`x` must be a numeric vector.")
    if (window %% 2 == 0 || window < 1)  stop("`window` must be a positive odd integer.")

    n <- length(x)
    half_w <- floor(window / 2)

    # Constant padding: repeat first and last values
    padded <- c(
      rep(x[1], half_w),
      x,
      rep(x[n], half_w)
    )

    # Allocate result
    result <- numeric(n)

    # Slide the window
    for (i in seq_len(n))
    {
      window_vals <- padded[(i-half_w):(i+half_w) + half_w]
      result[i] <- median(window_vals)
    }

    return(result)
  }

  for (id in unique(qsm$axis_ID))
  {
    axis = qsm[axis_ID == id]
    R = smooth_radius(axis$radius, fun = median)
    qsm$radius[qsm$axis_ID == id] = R
  }

  for (id in unique(qsm$axis_ID))
  {
    axis = qsm[axis_ID == id]
    R = smooth_radius(axis$radius, fun = mean)
    qsm$radius[qsm$axis_ID == id] = R
  }

  qsm = qsm_volume(qsm)

  return(qsm)
}


qsm_radius = function(qsm, tree, R0, tip_radius = 0.005, power = 1.1)
{
  plot = FALSE

  root = which(qsm$parent_ID == 0)
  w0 = qsm[["subtree_length"]][root]
  wi = qsm[["subtree_length"]]

  # Actual measurement for each cylindre using RANSAC if it work or falling back
  # to barycenter + average distance
  centroid = list()
  centroid$centroidX = (qsm$startX + qsm$endX)/2
  centroid$centroidY = (qsm$startY + qsm$endY)/2
  centroid$centroidZ = (qsm$startZ + qsm$endZ)/2
  centroid = as.data.frame(centroid)
  nearest = FNN::knnx.index(data = centroid, query = tree@data[,1:3], algorithm = "kd_tree", k = 1)

  xyz = tree@data[, 1:3]
  xyz$cyl_ID = qsm$cyl_ID[nearest]

  if (plot)
  {
    x = plot_qsm(qsm, cylinder = FALSE, add = c(0,0))
    cc = centroid
    rgl::points3d(cc, size = 3)
    rgl::points3d(xyz, col = lidR:::set.colors(xyz$cyl_ID, pal = lidR::random.colors(nrow(centroid))))
    lidR:::.pan3d(2)
  }

  data.table::setkey(qsm, cyl_ID)
  data.table::setkey(xyz, cyl_ID)

  r2 = rep(NA_real_, nrow(qsm))
  r_rescue = rep(NA_real_, nrow(qsm))

  fit_circle_least_squares <- function(xy)
  {
    x <- xy[,1]
    y <- xy[,2]
    A <- cbind(2*x, 2*y, rep(1, length(x)))
    b <- x^2 + y^2
    sol <- qr.solve(A, b)
    xc <- sol[1]
    yc <- sol[2]
    r  <- sqrt(sol[3] + xc^2 + yc^2)
    list(xc = xc, yc = yc, radius = r)
  }

  subs = split(xyz, by = "cyl_ID")
  axis = split(qsm, by = "cyl_ID")

  for (i in 1:nrow(qsm))
  {
    id = qsm$cyl_ID[i]

    if (id < 1) next

    id = as.character(id)

    axe = qsm[[id]]
    sub = subs[[id]]

    if (is.null(sub)) next

    sub = sub[,1:3]

    if (nrow(sub) < 100) next

    sub = as.matrix(sub)

    start = as.numeric(axe[,1:3])
    end = as.numeric(axe[,4:6])

    if (FALSE)
    {
      rgl::points3d(sub)
      rgl::segments3d(rbind(start, end), lwd = 3)
    }

    M = compute_rotation_matrix(start, end)
    sub = sub %*% t(M)

    if (nrow(sub) < 5) next

    r = lidRtls:::ransac_circle(sub, num_iterations = 100, inlier_threshold = 0.02)

    valid = is.valid.circle(r$radius, r$covered_arc_degree, r$percentage_inlier*100, r$percentage_inside*100)

    if (plot & valid)
    {
      col = if(valid) "darkgreen" else "red"
      title = paste0("ID = ",  i, " | arc = ", round(r$covered_arc_degree, 1),  " | inline = ", round(r$percentage_inlier*100), "% | R = ", round(r$radius, 2))
      plot(sub, asp = 1, pch = 19, cex = 0.5, main = title)
      symbols(r$center_x, r$center_y, circles = r$radius, add = TRUE, fg = col, inches = FALSE, lwd = 2)
      symbols(r$center_x, r$center_y, circles = r$radius+0.01, add = TRUE, fg = col, inches = FALSE, lty = 3, lwd = 2)
      symbols(r$center_x, r$center_y, circles = r$radius-0.01, add = TRUE, fg = col, inches = FALSE, lty = 3, lwd = 2)
    }

    if (!valid)
    {
      r_rescue[i] = fit_circle_least_squares(sub)$radius
    }
    else
    {
      #plot(sub, asp = 1, pch = 19, cex = 0.5)
      #symbols(r$center_x, r$center_y, circles = r$radius, add = TRUE, fg = "red", inches = FALSE)
      r2[i] = r$radius
    }

    #rgl::points3d(sub)
    #lidR::add_circle3d(c(0,0), xc, yc, r, mean(sub[,3]))
  }

  qsm[["radius"]] = r2
  qsm[["rescue_radius"]] = r_rescue
  qsm[["measure"]] = 0L
  qsm[["measure"]][!is.na(r2)] = 1L
  qsm$radius[qsm$subtree_length == 0] = 0

  if (plot)
  {
    x = plot_qsm(qsm, cylinder = TRUE)
    plot(tree, add = x, size = 2, pal = "brown")
  }

  # Maybe we have no bottom radius for this tree at all.
  # In this case we will use the rescue radius (least square method)
  # as replacement
  main_axis = qsm[axis_ID == 1]
  L = main_axis$subtree_length[1]
  first_ten_percent = main_axis[main_axis$subtree_length > L*0.90]
  if (all(is.na(first_ten_percent$radius)))
  {
    if (!all(is.na(first_ten_percent$rescue_radius)))
    {
      qsm$radius[qsm$cyl_ID %in% first_ten_percent$cyl_ID] = first_ten_percent$rescue_radius
      qsm$measure[qsm$cyl_ID %in% first_ten_percent$cyl_ID] = 2L
      first_ten_percent$radius = first_ten_percent$rescue_radius
    }
    else
    {

      idx = which(!is.na(main_axis$rescue_radius))

      if (length(idx) == 0)  stop("Internal error: no measurement on the bottom of the tree.")

      r0 = main_axis$rescue_radius[idx[1]]
      first_ten_percent$radius = r0
    }
  }
  qsm$rescue_radius = NULL

  if (plot)
  {
    plot_qsm(qsm, cylinder = TRUE)
    plot(tree, add = c(0,0), size = 2, pal = "brown")
  }

  # Interpolate NAs
  for (id in unique(qsm$axis_ID))
  {
    axis = qsm[axis_ID == id]
    idcyls = which(!is.na(axis$radius))
    idnacyls = which(is.na(axis$radius))

    if (length(idnacyls) == 0) next

    if (length(idcyls) > 0)
    {
      nearest = sapply(idnacyls, function(x) idcyls[which.min(abs(idcyls - x))])

      w0 = axis$subtree_length[nearest]
      wi = axis$subtree_length[idnacyls]
      radius = axis$radius[nearest]
      s = (wi / w0)^power
      r = s * radius + tip_radius
      if (any(is.infinite(r))) stop("Internal error: infinite radius detected")
      axis$radius[idnacyls] = r
      qsm$radius[qsm$axis_ID == id] = axis$radius
    }
  }

  if (plot)
  {
    x = plot_qsm(qsm, cylinder = TRUE)
    plot(tree, add = x, size = 2, pal = "brown")
  }

  R0 = sort(first_ten_percent$radius, TRUE)
  R0 = mean(R0[1:min(5, length(R0))])

  # AdTree allometric model
  qsm$radius_theorique = qsm_radius_adtree(qsm, tree, R0, tip_radius, power)$radius

  is_nas = is.nan(qsm$radius) | is.na(qsm$radius)
  qsm$radius[is_nas] = qsm$radius_theorique[is_nas]
  qsm$radius = ifelse(qsm$radius > 2*qsm$radius_theorique, qsm$radius_theorique, qsm$radius)
  qsm$radius_theorique = NULL

  # Smooth everything
  smooth_radius <- function(x, window = 3)
  {
    if (!is.numeric(x)) stop("x must be numeric")
    if (length(x) < 3) return(x)
    if (window %% 2 == 0) stop("window must be an odd number")

    half <- floor(window / 2)
    n <- length(x)

    # Pad the vector at both ends
    padded <- c(rep(x[1], half), x, rep(x[n], half))

    # Compute smoothed values
    result <- sapply(1:n, function(i) {
      median(padded[i:(i + 2 * half)])
    })

    return(result)
  }

  for (id in unique(qsm$axis_ID))
  {
    axis = qsm[axis_ID == id]
    R = smooth_radius(axis$radius, w = 5)
    qsm$radius[qsm$axis_ID == id] = R
  }

  qsm = qsm_volume(qsm)

  return(qsm)
}

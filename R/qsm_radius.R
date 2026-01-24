qsm_radius = function(qsm, tree, tip_radius = 0.0025)
{
  cat("Measuring diameters... ") ; ti = tic()

  H <- max(tree$hag)
  R0_first_guess <- DBH_vs_H_allometry(H)/2
  R0_second_guess <- find_root_radius(tree, qsm)
  R0 = max(R0_first_guess, R0_second_guess)

  qsm <- qsm_conic_allometry(qsm, 2*R0, tip_radius)
  conic <- qsm$radius

  if (R0 < 0.04) return(qsm)

  qsm <- qsm_measure(tree, qsm, sarc = 180, sins = 0.2, sinl = 0.3, srmeas = 0.03)
  #qsm <- qsm_measure_r(tree, qsm)
  qsm <- qsm_polynomial_fitting(qsm, tip_radius)

  # If we still have NAs on main axis it means interpolation failed on main
  # axis. We have no data.
  axis_ID <- NULL
  main_axis <- qsm[axis_ID == 1]
  if (anyNA(main_axis$radius))
  {
    warning("Not a single valid measure for this tree. The QSM is a pure reconstruction based on allometry", call. = FALSE)
    qsm <- qsm_conic_allometry(qsm, R0_first_guess, tip_radius)
    return(qsm)
  }

  qsm <- qsm_reconstruction_r(qsm, tip_radius)

  toc(ti)

  return(qsm)
}


qsm_conic_allometry = function(qsm, R0, tip_radius = 0.0025, power = 1)
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

qsm_measure = function(tree, qsm, sarc = 180, sins = 0.2, sinl = 0.3, srmeas = 0.03)
{
  # qsm_conic_allometry must be computed first
  qsm$theoric_radius = qsm$radius
  theoric_radius = qsm$radius
  qsm$radius = NULL

  qsm = qsm_measure_cpp(tree@data, qsm, sarc, sins, sinl, srmeas)
  data.table::setDT(qsm)
  return(qsm)
}

qsm_polynomial_fitting = function(qsm, tip_radius)
{
  qsm = qsm_polynomial_fitting_cpp(qsm, tip_radius)
  data.table::setDT(qsm)
  return(qsm[])
}

qsm_reconstruction = function(qsm, tip_radius)
{
  qsm = qsm_reconstruction_cpp(qsm, tip_radius)
  data.table::setDT(qsm)
  return(qsm[])
}

# *******************************************
# Original functions before to convert to C++
# *******************************************

qsm_measure_r = function(tree, qsm)
{
  # qsm_conic_allometry must be computed first
  qsm$theoric_radius = qsm$radius
  theoric_radius = qsm$radius
  qsm$radius = NULL

  cyl_ID <- NULL

  # Assign each point to an edge of the qsm graph
  centroid = list()
  centroid$centroidX = (qsm$startX + qsm$endX)/2
  centroid$centroidY = (qsm$startY + qsm$endY)/2
  centroid$centroidZ = (qsm$startZ + qsm$endZ)/2
  centroid = as.data.frame(centroid)
  nearest = FNN::knnx.index(data = centroid, query = tree@data[,1:3], algorithm = "kd_tree", k = 1)
  xyz = tree@data[, 1:3]
  xyz$cyl_ID = qsm$cyl_ID[nearest]

  if (FALSE)
  {
    x = plot_qsm(qsm, cylinder = FALSE, add = c(0,0))
    cc = centroid
    rgl::points3d(cc, size = 3)
    rgl::points3d(xyz, col = lidR:::set.colors(xyz$cyl_ID, pal = lidR::random.colors(nrow(centroid))))
    lidR:::.pan3d(2)
  }

  # Split data by cylinder (graph edges)
  data.table::setkey(qsm, cyl_ID)
  data.table::setkey(xyz, cyl_ID)

  r2 = rep(NA_real_, nrow(qsm))
  r_rescue = rep(NA_real_, nrow(qsm))

  subs = split(xyz, by = "cyl_ID",  keep.by = FALSE)
  starts = qsm[, c(1:3, 7)]
  ends = qsm[, c(4:6, 7)]
  starts = split(starts, by = "cyl_ID",  keep.by = FALSE)
  ends = split(ends, by = "cyl_ID",  keep.by = FALSE)
  theoric_radii = qsm[, c("cyl_ID", "theoric_radius")]
  theoric_radii = split(theoric_radii, by = "cyl_ID",  keep.by = FALSE)

  # Actual measurement for each cylindre using RANSAC if it works
  is.valid.circle = function(radius, angle_range, pinliner, pinside)
  {
    if (radius < 0.03)  return(FALSE)
    if (pinside > 20) return(FALSE)
    return(angle_range > 120 & pinliner > 30)
  }

  for (i in 1:nrow(qsm))
  {
    id = qsm$cyl_ID[i]

    if (id < 1) next # skip prolongations

    id = as.character(id)

    sub = subs[[id]]

    if (is.null(sub)) next

    theoric_radius = as.numeric(theoric_radii[[id]])

    # No need to try to measure anything if the theory is < 4 cm
    if (theoric_radius < 0.04) next

    sub = as.matrix(sub)

    start = as.numeric(starts[[id]])
    end = as.numeric(ends[[id]])

    if (FALSE)
    {
      rgl::points3d(sub)
      rgl::axes3d()
      rgl::segments3d(rbind(start, end), lwd = 3)
    }

    M = compute_rotation_matrix(start, end)
    sub = sub %*% t(M)

    if (nrow(sub) < 50) next

    r = ransac_circle(sub, num_iterations = 100, inlier_threshold = 0.02)

    valid = is.valid.circle(r$radius, r$covered_arc_degree, r$percentage_inlier*100, r$percentage_inside*100)

    if (valid)
    {
      r2[i] = r$radius

      if (FALSE)
      {
        col = if(valid) "darkgreen" else "red"
        title = paste0("ID = ",  i, " | arc = ", round(r$covered_arc_degree, 1),  " | inliner = ", round(r$percentage_inlier*100), "% insider = ", round(r$percentage_inside*100), "% | R = ", round(r$radius, 2))
        plot(sub, asp = 1, pch = 19, cex = 0.5, main = title)
        points(r$center_x, r$center_y, pch = 4)
        symbols(r$center_x, r$center_y, circles = r$radius, add = TRUE, fg = col, inches = FALSE, lwd = 2)
        symbols(r$center_x, r$center_y, circles = r$radius+0.01, add = TRUE, fg = col, inches = FALSE, lty = 3, lwd = 2)
        symbols(r$center_x, r$center_y, circles = r$radius-0.01, add = TRUE, fg = col, inches = FALSE, lty = 3, lwd = 2)
      }
    }
  }

  qsm[["radius"]] = data.table::copy(r2)

  if (FALSE)
  {
    x = plot_qsm(qsm, cylinder = TRUE)
    plot(tree, add = x, size = 2, pal = "brown")
  }

  return(qsm)
}

qsm_polynomial_fitting_r = function(qsm, tip_radius)
{
  radius <- axis_ID <- NULL

  # We now have some cynlindred measure with ransac
  # Polynomial interpolation for each axis with enought measurements
  for (i in sort(unique(qsm$axis_ID)))
  {
    axe = qsm[axis_ID == i]
    if (sum(!is.na(axe$radius)) <= 6) next

    mes = axe[!is.na(radius)]
    mod <- stats::nls(
      radius ~ tip_radius + a * subtree_length + b * subtree_length^2,
      data = mes,
      start = list(a = 0.01, b = 0.001),
      algorithm = "port",
      lower = c(-Inf, -Inf),
      upper = c(Inf, Inf)   # b > 0
    )

    r = stats::predict(mod, axe)
    qsm[axis_ID == i, radius := r]

    if (FALSE)
    {
      # plot(axe$theoric_radius, max(axe$subtree_length) -axe$subtree_length,
      #      asp = 0.01,
      #      main = NULL,
      #      pch  = 19,
      #      cex = 0.5,
      #      col = "gray",
      #      xlab = "Radius (m)",
      #      ylab = "Distance to root (m)",
      #      xlim = c(0, 0.2))
      plot(axe$radius, max(axe$subtree_length) -axe$subtree_length,
           pch = 19, cex = 0.6,
           xlab = "Radius (m)",
           ylab = "Distance to root (m)",
           xlim = c(0, max(r, na.omit(axe$radius))))
      points(r, max(axe$subtree_length) -axe$subtree_length, col = "blue", pch = 19, cex = 0.5)

      legend("topright",
             legend = c("Theoric radius", "Measured radius", "Interpolation"),
             col = c("gray", "black", "blue"),
             pch = 19,
             pt.cex = c(0.5, 0.6, 0.5))
    }
  }

  return(qsm[])
}

qsm_reconstruction_r = function(qsm, tip_radius)
{
  axis_ID <- radius <- branch_order <- NULL

  # We have some cylinder measured and interpolated. A large portion of the tree is missing.
  # We loop on each axes by branch order. If we have some NAs we compare to the theory.
  # We ensure the theory is not bigger than the parent branch. If the theory is bigger
  # We rescale to be the size of the parent
  # And if we have enough measurement we scale the theory to something realistic

  # Loop by branch order
  for (i in sort(unique(qsm$branch_order)))
  {
    if (i == 1) next

    order = qsm[branch_order == i]
    axis_IDs = unique(order$axis_ID)

    for (j in axis_IDs)
    {
      axe = order[axis_ID == j]

      # Find parent radius and subtree length
      parent = axe$parent_ID[1]
      parent_id = which(qsm$cyl_ID == parent)
      if (length(parent_id) == 0) stop("Internal error. No parent.")

      # Compute the theoretical radius
      r0 = qsm$radius[parent_id]
      w0 = qsm$subtree_length[parent_id]
      wi = axe$subtree_length
      s  = (wi / w0)^1.1
      r  = tip_radius + s * (r0 - tip_radius)

      # We have only NAs. We can only use the theory
      if (all(is.na(axe$radius)))
      {

      }
      # We have some measures. Use a comparison between measure and theory to rescale theory
      else if(anyNA(axe$radius))
      {
        r_mes = axe$radius

        # We need at least 3 measurements otherwise it is easy to have
        # false ransac fitting. Other wise use radius theoretic
        if (sum(!is.na(r_mes)) >= 3)
        {
          # Look at the average ratio between the theory in and the measures.
          # If measures tends to be bigger increase theory. Or opposite
          ratio = stats::median(r_mes/r, na.rm = TRUE)
          r = tip_radius + s * (r0*ratio - tip_radius)
        }
      }

      qsm[branch_order == i & axis_ID == j, radius := r]
    }
  }

  if (FALSE)
  {
    x = plot_qsm(qsm, cylinder = TRUE)
    plot(tree, add =x, size = 2, pal = "brown")
  }

  return(qsm[])
}

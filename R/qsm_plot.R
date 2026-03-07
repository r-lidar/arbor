#' Plot QSM Data in 3D
#'
#' Plots the QSM data as a series of connected cylinders, with color-coding based on branch attributes.
#'
#' @param x A QSM of QSF
#' @param qsm A QSM
#' @param add A numeric vector for translation offsets. Like in the lidR package.
#' @param sides Number of sides for each cylinder.
#' @param color The attribute for color mapping.
#' @param pal Color palette
#' @param ... Unused (for S3 compatibility).
#' @param skeleton,cylinder boolean. plot the skeleton or the cylinders or both.
#' @method plot qsm
#' @export
#' @rdname plot
#' @md
plot.qsm = function(x, ...)
{
  plot_qsm(x, ...)
}

#' @export
#' @rdname plot
plot_qsm = function(qsm, add = NULL, sides = 12, color = "cyl_ID", skeleton = TRUE, cylinder = TRUE, pal = c("blue", "green", "yellow", "orange", "red"), ...)
{
  # --- Input Validation ---
  if (!is.data.frame(qsm)) {
    stop("Input must be a data.frame")
  }

  color_palette <- grDevices::colorRampPalette(pal)

  # --- Global Translation Logic ---
  # translateXYZ is for compatibility with computree
  if (!is.null(add)) {
    tx = add[1]; ty = add[2]; tz = 0
  } else {
    if (!"translateX" %in% names(qsm)) tx = min(qsm$startX) else tx = qsm$translateX[1]
    if (!"translateY" %in% names(qsm)) ty = min(qsm$startY) else ty = qsm$translateY[1]
    if (!"translateZ" %in% names(qsm)) tz = 0 else tz = qsm$translateZ[1]

    rgl::open3d()
  }

  # --- Check Cylinder Capability ---
  local_cylinder <- cylinder
  if (!"radius" %in% names(qsm)) local_cylinder <- FALSE

  # --- Color Logic ---
  colattr = qsm[[color]]
  if (color %in% names(qsm)) {
    if (color == "branch_order") {
      colattr = colattr
    } else if (is.logical(colattr)) {
      colattr = colattr + 1
    } else {
      # Handle case where min == max to avoid seq error
      mn <- min(colattr, na.rm = TRUE)
      mx <- max(colattr, na.rm = TRUE)
      if (mn == mx) {
        colattr <- rep(1, length(colattr))
      } else {
        colattr = findInterval(colattr, seq(mn, mx, length.out = 20))
      }
    }
    colattr[colattr < 1] <- 1
    colors_mapped <- color_palette(max(colattr, na.rm=TRUE))[colattr]
  } else {
    colors_mapped = rep("black", nrow(qsm))
  }

  # --- Apply Translation ---
  qsm$startX <- qsm$startX - tx
  qsm$startY <- qsm$startY - ty
  qsm$startZ <- qsm$startZ - tz
  qsm$endX   <- qsm$endX - tx
  qsm$endY   <- qsm$endY - ty
  qsm$endZ   <- qsm$endZ - tz

  # Update bounds for axes
  z_min_glob <- min(qsm$startZ, qsm$endZ)
  z_max_glob <- max(qsm$startZ, qsm$endZ)

  # --- Generate and Render Mesh ---
  if (local_cylinder) {
    mesh <- cylinders_as_mesh(qsm, sides = sides, color_vec = colors_mapped)
    rgl::shade3d(mesh)
  }

  # --- Render Skeleton ---
  if (skeleton) {
    # Interleave Start and End for segments
    pts <- matrix(NA, nrow=nrow(qsm)*2, ncol=3)
    pts[seq(1, nrow(pts), 2), ] <- as.matrix(qsm[, c("startX","startY","startZ")])
    pts[seq(2, nrow(pts), 2), ] <- as.matrix(qsm[, c("endX","endY","endZ")])

    line_cols <- rep(colors_mapped, each=2)

    rgl::segments3d(pts, col = line_cols)
    rgl::points3d(as.matrix(qsm[, c("startX","startY","startZ")]), col = colors_mapped)
  }

  # --- Axes ---
  z_ticks <- seq(floor(z_min_glob), ceiling(z_max_glob), by = 0.5)
  rgl::axis3d("z", at = z_ticks, labels = as.character(z_ticks), col = "black")
  rgl::axis3d("x", col = "black")
  rgl::axis3d("y", col = "black")

  lidR:::.pan3d(2)

  return(invisible(c(tx, ty)))
}

cylinders_as_mesh <- function(qsm, sides = 16, color_vec = "black")
{
  mesh_data <- qsm_mesh_cpp(qsm, sides)

  n_cyl <- nrow(qsm)
  n_verts_total <- ncol(mesh_data$vertices)

  id = match(mesh_data$NodeID, qsm$cyl_ID)
  Final_Colors = color_vec[id]

  mesh <- rgl::qmesh3d(
    vertices = mesh_data$vertices,
    indices  = mesh_data$indices,
    homogeneous = FALSE
  )

  mesh$material$color <- Final_Colors

  return(mesh)
}


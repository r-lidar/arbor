#' Plot QSM Data in 3D
#'
#' Plots the QSM data as a series of connected cylinders, with color-coding based on branch attributes.
#'
#' @param qsm A data frame containing QSM data with segment attributes.
#' @param add A numeric vector for translation offsets. Like in the lidR package.
#' @param sides Number of sides for each cylinder.
#' @param color The attribute for color mapping.
#' @param ... Unused (for S3 compatibility).
#' @param skeleton,cylinder boolean. plot the skeleton or the cylinders or both.
#' @method plot qsm
#' @export
#' @md
plot.qsm = function(qsm, add = NULL, sides = 12, color = "cyl_ID", skeleton = TRUE, cylinder = TRUE, pal = c("blue", "green", "yellow", "orange", "red"), ...)
{
  color_palette <- grDevices::colorRampPalette(pal)
  colattr = qsm[[color]]

  if (!"radius" %in% names(qsm)) cylinder = FALSE

  # --- Color Logic (Kept same) ---
  if (color %in% names(qsm)) {
    if (color == "branch_order") {
      colattr = colattr
    } else if (is.logical(colattr)) {
      colattr = colattr + 1
    } else {
      colattr = findInterval(colattr, seq(min(colattr), max(colattr), length.out = 20))
    }
    # Ensure indices are within bounds of palette
    colattr[colattr < 1] <- 1
    colors_mapped <- color_palette(max(colattr, na.rm=TRUE))[colattr]
  } else {
    colors_mapped = rep("black", nrow(qsm))
  }

  # --- Translation Logic ---
  if (!"translateX" %in% names(qsm)) qsm$translateX = min(qsm$startX)
  if (!"translateY" %in% names(qsm)) qsm$translateY = min(qsm$startY)
  if (!"translateZ" %in% names(qsm)) qsm$translateZ = 0

  if (!is.null(add)) {
    tx = add[1]; ty = add[2]; tz = 0
  } else {
    tx = qsm$translateX[1]
    ty = qsm$translateY[1]
    tz = 0
    rgl::open3d()
  }

  # Apply Translation to the dataframe temporarily for plotting
  qsm_t <- qsm
  qsm_t$startX <- qsm$startX - tx
  qsm_t$startY <- qsm$startY - ty
  qsm_t$startZ <- qsm$startZ - tz
  qsm_t$endX   <- qsm$endX - tx
  qsm_t$endY   <- qsm$endY - ty
  qsm_t$endZ   <- qsm$endZ - tz

  # --- FAST RENDERING ---
  if (cylinder)
  {
    # Generate one single mesh for all cylinders
    mesh <- cylinders_as_mesh(qsm_t, sides = sides, color_vec = colors_mapped)
    rgl::shade3d(mesh)
  }

  # --- Skeleton (Lines) ---
  if (skeleton)
  {
    # Fast segment plotting
    # Interleave Start and End for segments3d
    pts <- matrix(NA, nrow=nrow(qsm)*2, ncol=3)
    pts[seq(1, nrow(pts), 2), ] <- as.matrix(qsm_t[, c("startX","startY","startZ")])
    pts[seq(2, nrow(pts), 2), ] <- as.matrix(qsm_t[, c("endX","endY","endZ")])

    # Interleave colors
    line_cols <- rep(colors_mapped, each=2)

    rgl::segments3d(pts, col = line_cols)

    # Points at start
    rgl::points3d(qsm_t[, c("startX","startY","startZ")], col = colors_mapped)
  }

  # --- Axes ---
  if (is.null(add))
  {
    z_limits <- range(c(qsm_t$startZ, qsm_t$endZ))
    z_ticks <- seq(floor(z_limits[1]), ceiling(z_limits[2]), by = 0.5)
    rgl::axis3d("z", at = z_ticks, labels = as.character(z_ticks), col = "black")
    rgl::axis3d("x", col = "black")
    rgl::axis3d("y", col = "black")
  }

  if (is.null(add))
  {
    lidR:::.pan3d(2)
  }

  return(invisible(c(tx, ty)))
}

cylinders_as_mesh <- function(qsm, sides = 16, color_vec = "black")
{
  n_cyl <- nrow(qsm)

  # 1. Prepare Data Matrices
  P1 <- as.matrix(qsm[, c("startX", "startY", "startZ")])
  P2 <- as.matrix(qsm[, c("endX", "endY", "endZ")])
  R  <- qsm$radius

  # 2. Compute Cylinder Axes
  Vz <- P2 - P1
  h  <- sqrt(rowSums(Vz^2))
  Vz <- Vz / h

  # 3. Compute Local Coordinate Systems
  Vtemp <- matrix(c(0,0,1), nrow = n_cyl, ncol = 3, byrow = TRUE)

  # Check for parallel vectors
  is_parallel <- abs(Vz[,3]) > 0.99

  # FIX: Correctly assign (1,0,0) to all matching rows
  if (any(is_parallel)) {
    Vtemp[is_parallel, 1] <- 1
    Vtemp[is_parallel, 2] <- 0
    Vtemp[is_parallel, 3] <- 0
  }

  # Cross product to get local X axis
  Vx <- cbind(Vz[,2]*Vtemp[,3] - Vz[,3]*Vtemp[,2],
              Vz[,3]*Vtemp[,1] - Vz[,1]*Vtemp[,3],
              Vz[,1]*Vtemp[,2] - Vz[,2]*Vtemp[,1])

  Vx <- Vx / sqrt(rowSums(Vx^2))

  # Cross product to get local Y axis
  Vy <- cbind(Vz[,2]*Vx[,3] - Vz[,3]*Vx[,2],
              Vz[,3]*Vx[,1] - Vz[,1]*Vx[,3],
              Vz[,1]*Vx[,2] - Vz[,2]*Vx[,1])

  # 4. Generate Vertices
  theta <- seq(0, 2*pi, length.out = sides + 1)[-(sides + 1)]
  cos_t <- cos(theta)
  sin_t <- sin(theta)

  R_long  <- rep(R, each = sides)
  Vx_long <- Vx[rep(1:n_cyl, each = sides), ]
  Vy_long <- Vy[rep(1:n_cyl, each = sides), ]

  cos_long <- rep(cos_t, times = n_cyl)
  sin_long <- rep(sin_t, times = n_cyl)

  Offset <- R_long * (Vx_long * cos_long + Vy_long * sin_long)

  P1_long <- P1[rep(1:n_cyl, each = sides), ]
  P2_long <- P2[rep(1:n_cyl, each = sides), ]

  Verts_bot <- P1_long + Offset
  Verts_top <- P2_long + Offset
  Vertices  <- rbind(Verts_bot, Verts_top)

  # 5. Generate Indices
  i <- 1:sides
  j <- c(2:sides, 1)

  total_verts_bot <- n_cyl * sides
  offsets <- (0:(n_cyl-1)) * sides

  idx_row1 <- rep(i, n_cyl) + rep(offsets, each=sides)
  idx_row2 <- rep(j, n_cyl) + rep(offsets, each=sides)
  idx_row3 <- rep(j, n_cyl) + rep(offsets, each=sides) + total_verts_bot
  idx_row4 <- rep(i, n_cyl) + rep(offsets, each=sides) + total_verts_bot

  Indices <- rbind(idx_row1, idx_row2, idx_row3, idx_row4)

  # 6. Handle Colors
  if (length(color_vec) == n_cyl) {
    col_expanded <- rep(color_vec, each = sides)
    Final_Colors <- c(col_expanded, col_expanded)
  } else {
    Final_Colors <- rep(color_vec[1], nrow(Vertices))
  }

  # 7. Create mesh3d Object
  mesh <- rgl::qmesh3d(
    vertices = t(Vertices),
    indices  = Indices,
    homogeneous = FALSE
  )

  mesh$material$color <- Final_Colors

  return(mesh)
}


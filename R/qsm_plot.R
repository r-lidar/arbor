#' Plot QSM Data in 3D
#'
#' Plots the QSM data as a series of connected cylinders, with color-coding based on branch attributes.
#'
#' @param x,qsm A data frame containing QSM data with segment attributes.
#' @param add A numeric vector for translation offsets. Like in the lidR package.
#' @param sides Number of sides for each cylinder.
#' @param color The attribute for color mapping.
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
  # --- Input Normalization ---
  # Ensure qsm is a list of dataframes
  if (is.data.frame(qsm)) {
    qsm_list <- list(qsm)
  } else if (is.list(qsm)) {
    qsm_list <- qsm
  } else {
    stop("Input must be a data.frame or a list of data.frames")
  }

  color_palette <- grDevices::colorRampPalette(pal)

  # --- Global Translation Logic ---
  # Determine translation based on 'add' or the first QSM object
  if (!is.null(add)) {
    tx = add[1]; ty = add[2]; tz = 0
  } else {
    # Check first QSM for defaults
    q1 <- qsm_list[[1]]
    if (!"translateX" %in% names(q1)) tx = min(q1$startX) else tx = q1$translateX[1]
    if (!"translateY" %in% names(q1)) ty = min(q1$startY) else ty = q1$translateY[1]
    if (!"translateZ" %in% names(q1)) tz = 0 else tz = q1$translateZ[1]

    rgl::open3d()
  }

  # Storage for batch rendering
  meshes_to_merge   <- list()
  skeleton_segments <- list()
  skeleton_colors   <- list()
  skeleton_points   <- list()
  skeleton_pcolors  <- list()

  # Track global Z range for axes
  z_min_glob <- Inf
  z_max_glob <- -Inf

  # --- Loop over Objects ---
  for (i in seq_along(qsm_list))
  {
    q <- qsm_list[[i]]

    # Check per-object cylinder capability
    local_cylinder <- cylinder
    if (!"radius" %in% names(q)) local_cylinder <- FALSE

    # --- Color Logic ---
    colattr = q[[color]]
    if (color %in% names(q)) {
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
      colors_mapped = rep("black", nrow(q))
    }

    # --- Apply Translation ---
    q_t <- q
    q_t$startX <- q$startX - tx
    q_t$startY <- q$startY - ty
    q_t$startZ <- q$startZ - tz
    q_t$endX   <- q$endX - tx
    q_t$endY   <- q$endY - ty
    q_t$endZ   <- q$endZ - tz

    # Update bounds for axes
    if (is.null(add)) {
      z_min_glob <- min(z_min_glob, min(q_t$startZ), min(q_t$endZ))
      z_max_glob <- max(z_max_glob, max(q_t$startZ), max(q_t$endZ))
    }

    # --- Generate Mesh ---
    if (local_cylinder) {
      # Assuming cylinders_as_mesh returns a mesh3d object
      mesh <- cylinders_as_mesh(q_t, sides = sides, color_vec = colors_mapped)
      meshes_to_merge[[length(meshes_to_merge) + 1]] <- mesh
    }

    # --- Prepare Skeleton ---
    if (skeleton) {
      # Interleave Start and End for segments
      pts <- matrix(NA, nrow=nrow(q)*2, ncol=3)
      pts[seq(1, nrow(pts), 2), ] <- as.matrix(q_t[, c("startX","startY","startZ")])
      pts[seq(2, nrow(pts), 2), ] <- as.matrix(q_t[, c("endX","endY","endZ")])

      line_cols <- rep(colors_mapped, each=2)

      skeleton_segments[[length(skeleton_segments)+1]] <- pts
      skeleton_colors[[length(skeleton_colors)+1]]     <- line_cols

      skeleton_points[[length(skeleton_points)+1]]     <- as.matrix(q_t[, c("startX","startY","startZ")])
      skeleton_pcolors[[length(skeleton_pcolors)+1]]   <- colors_mapped
    }
  }

  # --- Fast Rendering ---

  # 1. Plot merged mesh
  if (length(meshes_to_merge) > 0) {
    if (length(meshes_to_merge) == 1) {
      rgl::shade3d(meshes_to_merge[[1]])
    } else {
      # Merge all meshes into a single mesh object for efficiency
      full_mesh <- do.call(merge, meshes_to_merge)
      rgl::shade3d(full_mesh)
    }
  }

  # 2. Plot skeleton (batched)
  if (skeleton && length(skeleton_segments) > 0) {
    full_segments <- do.call(rbind, skeleton_segments)
    full_seg_cols <- unlist(skeleton_colors)
    rgl::segments3d(full_segments, col = full_seg_cols)

    full_pts <- do.call(rbind, skeleton_points)
    full_pts_cols <- unlist(skeleton_pcolors)
    rgl::points3d(full_pts, col = full_pts_cols)
  }

  # --- Axes ---
  # Handle case with no data gracefully
  if (is.infinite(z_min_glob)) { z_min_glob <- 0; z_max_glob <- 1 }

  z_ticks <- seq(floor(z_min_glob), ceiling(z_max_glob), by = 0.5)
  rgl::axis3d("z", at = z_ticks, labels = as.character(z_ticks), col = "black")
  rgl::axis3d("x", col = "black")
  rgl::axis3d("y", col = "black")

  lidR:::.pan3d(2)

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


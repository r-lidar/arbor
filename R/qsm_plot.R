# @file qsm_plot.R
# Project: Arbor
#
# Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' Plot QSM Data in 3D
#'
#' Plots the QSM data as a series of connected cylinders, with color-coding based on branch attributes.
#'
#' @param x A QSM of QSF
#' @param qsm A QSM
#' @param add A numeric vector for translation offsets. Like in the lidR package.
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


# Helper function to ensure identical color mapping in both functions
map_qsm_colors <- function(colattr, pal, is_categorical) {
  if (is_categorical) {
    # Categorical logic: Direct indexing with modulo wrap-around
    colattr[is.na(colattr)] <- 1
    colattr[colattr < 1] <- 1
    idx <- pmin(as.integer(colattr), length(pal))
    return(pal[idx])
  } else {
    # Continuous logic: Fixed-bin gradient mapping
    mn <- min(colattr, na.rm = TRUE)
    mx <- max(colattr, na.rm = TRUE)

    if (mn == mx) {
      bins <- rep(1, length(colattr))
      n_bins <- 1
    } else {
      n_bins <- 100 # Use a fixed resolution for the ramp
      bins <- findInterval(colattr, seq(mn, mx, length.out = n_bins))
    }

    color_palette <- grDevices::colorRampPalette(pal)
    return(color_palette(n_bins)[bins])
  }
}

#' @export
#' @rdname plot
plot_qsm = function(qsm, add = NULL, color = "branch_order", skeleton = TRUE, cylinder = TRUE, pal = "auto", ...) {
  default_pal = c("blue", "green", "yellow", "orange", "red")
  branch_order_pal <- c("#552203","#ad3a01","#e59b16","#fdd63b","#9acd56", "#59ad4b","#1e490e")
  quality_pal = c("blue", "green", "yellow", "orange", "red")

  if (!is.data.frame(qsm)) stop("Input must be a data.frame")

  if (identical(pal, "auto")) {
    pal <- switch(color, branch_order = branch_order_pal, quality = quality_pal, default_pal)
  }

  if (!is.null(add)) {
    tx = add[1]; ty = add[2]; tz = 0
  } else {
    tx = min(qsm$startX); ty = min(qsm$startY); tz = 0
    rgl::open3d()
  }

  # Apply Translation
  qsm$startX <- qsm$startX - tx; qsm$startY <- qsm$startY - ty; qsm$startZ <- qsm$startZ - tz
  qsm$endX   <- qsm$endX - tx;   qsm$endY   <- qsm$endY - ty;   qsm$endZ   <- qsm$endZ - tz

  # Generate Colors (Centralized Logic)
  if (color %in% names(qsm)) {
    is_cat <- color %in% c("branch_order", "quality")
    colors_mapped <- map_qsm_colors(qsm[[color]], pal, is_cat)
  } else {
    colors_mapped <- rep("black", nrow(qsm))
  }

  # Render Mesh
  if (cylinder && "radius" %in% names(qsm)) {
    # Pass the already computed colors to avoid re-calculation mismatch
    mesh <- as_mesh(qsm, color, pal, precomputed_colors = colors_mapped)
    rgl::shade3d(mesh)
  }

  # Render Skeleton
  if (skeleton && nrow(qsm) > 0) {
    pts <- matrix(NA, nrow = nrow(qsm) * 2, ncol = 3)
    pts[seq(1, nrow(pts), 2), ] <- as.matrix(qsm[, c("startX", "startY", "startZ")])
    pts[seq(2, nrow(pts), 2), ] <- as.matrix(qsm[, c("endX", "endY", "endZ")])

    rgl::segments3d(pts, col = rep(colors_mapped, each = 2))
    rgl::points3d(as.matrix(qsm[, c("startX", "startY", "startZ")]), col = colors_mapped)
  }

  lidR:::.pan3d(2)
  return(invisible(c(tx, ty)))
}

as_mesh <- function(qsm, color = "cyl_ID", pal = c("blue", "green", "yellow", "orange", "red"), precomputed_colors = NULL) {
  if (nrow(qsm) == 0) return(NULL)

  # Use precomputed colors if provided, otherwise compute using the same logic
  if (!is.null(precomputed_colors)) {
    colors_mapped <- precomputed_colors
  } else {
    if (color %in% names(qsm)) {
      is_cat <- color %in% c("branch_order", "quality")
      colors_mapped <- map_qsm_colors(qsm[[color]], pal, is_cat)
    } else {
      colors_mapped <- rep("black", nrow(qsm))
    }
  }

  mesh_data <- qsm_mesh_cpp(qsm, 16)

  # Align colors with mesh vertices using NodeID (CylID)
  id = match(mesh_data$NodeID, qsm$cyl_ID)
  Final_Colors = colors_mapped[id]

  mesh <- rgl::qmesh3d(
    vertices = mesh_data$vertices,
    indices  = mesh_data$indices,
    homogeneous = FALSE
  )

  mesh$material$color <- Final_Colors
  mesh$material$specular <- "black"
  mesh$material$shininess <- 0

  return(mesh)
}

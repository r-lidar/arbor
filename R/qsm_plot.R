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
plot_qsm <- function(x, add = NULL, color = "branch_order", skeleton = TRUE, cylinder = TRUE, pal = "auto", ...) {
  if (!is.data.frame(x)) stop("Input must be a data.frame")

  offset <- add
  if (is.null(offset)) {
    offset <- c(min(x$startX), min(x$startY))
    rgl::open3d()
  }

  pal <- get_qsm_pal(color, pal)

  # Generate Colors for both Mesh and Skeleton
  if (color %in% names(x)) {
    is_cat <- color %in% c("branch_order", "quality")
    colors_mapped <- map_qsm_colors(x[[color]], pal, is_cat)
  } else {
    colors_mapped <- rep("black", nrow(x))
  }

  # Render Mesh
  if (cylinder && "radius" %in% names(x)) {
    mesh <- as_mesh.qsm(x, offset = offset, color = color, pal = pal, precomputed_colors = colors_mapped)
    rgl::shade3d(mesh)
  }

  # Render Skeleton
  if (skeleton && nrow(x) > 0) {
    # We shift coordinates here strictly for the skeleton plot, since 
    # as_mesh handles its own shifting internally for the mesh.
    pts <- matrix(NA, nrow = nrow(x) * 2, ncol = 3)
    pts[seq(1, nrow(pts), 2), 1] <- x$startX - offset[1]
    pts[seq(1, nrow(pts), 2), 2] <- x$startY - offset[2]
    pts[seq(1, nrow(pts), 2), 3] <- x$startZ
    
    pts[seq(2, nrow(pts), 2), 1] <- x$endX - offset[1]
    pts[seq(2, nrow(pts), 2), 2] <- x$endY - offset[2]
    pts[seq(2, nrow(pts), 2), 3] <- x$endZ

    rgl::segments3d(pts, col = rep(colors_mapped, each = 2))
    
    start_pts <- cbind(x$startX - offset[1], x$startY - offset[2], x$startZ)
    rgl::points3d(start_pts, col = colors_mapped)
  }

  lidR:::.pan3d(2)
  return(invisible(offset))
}
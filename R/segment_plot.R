# @file segment_plot.R
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

#' Plot functions
#'
#' `colorize_trees` assigns a RGB value per tree for rendering in third
#' party softwares.
#'
#' @param las LAS object
#' @param dtm SpatRaster DTM
#' @param ... propagated to lidR::plot
#' @param th threshold to remove points with few passages
#' @export
#' @rdname plot
plot_likelihood = function(las, dtm = NULL, ...)
{
  x <- lidR::plot(las, color = "pwood", legend = T, breaks = "quantile", ...)
  if (!is.null(dtm)) lidR::add_dtm3d(x, dtm)
  return(invisible(x))
}

#' @export
#' @rdname plot
plot_semantic = function(las, dtm = NULL, ...)
{
  x <- lidR::plot(las, color = "foliage", pal = foliage.colors, ...)
  if (!is.null(dtm)) lidR::add_dtm3d(x, dtm)
  return(invisible(x))
}

#' @export
#' @rdname plot
plot_instance = function(las, dtm = NULL, ...)
{
  set.seed(42)
  if ("UserData" %in% names(las)) las@data$treeID[las$UserData != ARBORTREE] = NA_integer_
  las@data$treeID[las@data$treeID == 0] = NA_integer_
  x <- lidR::plot(las, color = "treeID", ...)
  if (!is.null(dtm)) lidR::add_dtm3d(x, dtm)
  return(invisible(x))
}

#' @export
#' @rdname plot
plot_semantic_instance = function(las, dtm = NULL, ...)
{
  p1 <- las@data$foliage    # 0 = wood, 1/2 = foliage
  p2 <- las@data$treeID     # tree IDs (non-continuous)
  p2[p2 <= 0] <- NA_integer_   # 0 = no tree. A 0 index would drop rows in pal[p2, ]
  n <- max(p2, na.rm = TRUE)   # number of unique trees

  # Example pastel palette
  pal <- lidR::pastel.colors(n)   # or any vector of colors per tree
  pal <- t(grDevices::col2rgb(pal))

  # Create vector to store RGB for each point
  cols <- pal[p2, ]

  # Darken foliage points
  foliage = p1 >= 1
  cols[foliage,] =   cols[foliage,] * 0.7
  R = as.integer(cols[,1])
  G = as.integer(cols[,2])
  B = as.integer(cols[,3])
  R[is.na(R)] = 150L
  G[is.na(G)] = 150L
  B[is.na(B)] = 150L

  # Assign to LAS
  las@data$R = R
  las@data$G = G
  las@data$B = B

  free(cols)

  x <- lidR::plot(las, color = "RGB")
  if (!is.null(dtm)) lidR::add_dtm3d(x, dtm)
  return(invisible(x))
}

#' @export
#' @rdname plot
plot_passage = function(las, dtm = NULL, th = 0, ...)
{
  passage <- NULL
  passage <- lidR::filter_poi(las, passage > th)
  passage@data$passage <- log(passage$passage)
  x <- lidR::plot(passage, color = "passage", legend = T, ...)
  if (!is.null(dtm)) lidR::add_dtm3d(x, dtm)
  return(invisible(x))
}

#' @export
#' @rdname plot
#' @param darken_foliage boolean. In addition to one color per tree, the foliage is darkened
#' the same color that its corresponding tree to display semantic and instance segmentation in
#' one rendering.
colorize_trees = function(las, darken_foliage = TRUE)
{
  # Memory: avoid lidR::add_lasrgb() and las@data$x <- y. Both copy the whole data.table.
  # Columns are shared by reference and only R, G, B (and a negated treeID) are allocated.
  n = lidR::npoints(las)
  grey = 150L * 257L

  tid = las$treeID
  if ("UserData" %in% names(las))
  {
    # negative values for treeID such as C++ colorization skips them. Copy on modify: las is untouched.
    skip = which(las$UserData != ARBORTREE)
    tid[skip] = -tid[skip]
  }

  R = rep(grey, n)
  G = rep(grey, n)
  B = rep(grey, n)

  # Lightweight data.frame view on existing columns for C++ (no copy). C++ writes R, G, B in place.
  core = list(X = las$X, Y = las$Y, Z = las$Z, treeID = tid, R = R, G = G, B = B)
  if ("foliage" %in% names(las)) core$foliage = las$foliage
  data.table::setDT(core)
  colorize_trees_cpp(core, darken_foliage)

  # Shallow copy of the table so that adding columns by reference does not alter the input las
  data = data.table::setDT(as.list(las@data))
  data.table::set(data, j = c("R", "G", "B"), value = list(R, G, B))
  las@data = data

  format = las@header@PHB[["Point Data Format ID"]]
  if (!format %in% c(2L, 3L, 8L))
  {
    if ("NIR" %in% names(las)) format = 8L
    else if ("gpstime" %in% names(las)) format = 3L
    else format = 2L
    las@header@PHB[["Point Data Format ID"]] = format
  }

  return(las)
}





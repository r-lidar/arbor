#' Plot function
#'
#' Various useful plot functions. Colorize tree does not plot be assigns a RGB value per tree
#' for rendering.
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
  las@data$treeID[las$treeID < 0] = NA_integer_
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
  R = rep(150L, lidR::npoints(las))
  G = rep(150L, lidR::npoints(las))
  B = rep(150L, lidR::npoints(las))
  las = lidR::add_lasrgb(las, R, G, B)
  colorize_trees_cpp(las@data, darken_foliage)
  return(las)
}





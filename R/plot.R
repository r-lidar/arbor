#' @export
plot_anisotropy = function(las, dtm = NULL, ...)
{
  x <- lidR::plot(las, color = "anisotropy", legend = T, breaks = "quantile", ...)
  if (!is.null(dtm)) lidR::add_dtm3d(x, dtm)
  return(invisible(x))
}

#' @export
plot_semantic = function(las, dtm = NULL, ...)
{
  x <- lidR::plot(las, color = "foliage", pal = foliage.colors, ...)
  if (!is.null(dtm)) lidR::add_dtm3d(x, dtm)
  return(invisible(x))
}

#' @export
plot_instance = function(las, dtm = NULL, ...)
{
  x <- lidR::plot(las, color = "treeID", ...)
  if (!is.null(dtm)) lidR::add_dtm3d(x, dtm)
  return(invisible(x))
}

#' @export
plot_semantic_instance = function(las, dtm = NULL, ...)
{
  p1 <- las@data$foliage    # 0 = wood, 1/2 = foliage
  p2 <- las@data$treeID     # tree IDs (non-continuous)
  n <- max(p2, na.rm = TRUE)   # number of unique trees

  # Example pastel palette
  pal <- pastel.colors(n)   # or any vector of colors per tree
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
plot_passage = function(las, dtm = NULL)
{
  passage <- NULL
  passage <- lidR::filter_poi(las, passage > 1)
  passage@data$passage <- log(passage$passage)
  x <- lidR::plot(passage, color = "passage", legend = T)
  if (!is.null(dtm)) lidR::add_dtm3d(x, dtm)
  return(invisible(x))
}





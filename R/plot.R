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
  n <- length(unique(p2))   # number of unique trees

  # Example pastel palette
  pal <- pastel.colors(n)   # or any vector of colors per tree

  # Map tree IDs to palette colors
  id_unique <- sort(unique(p2))
  id_to_color <- setNames(pal, id_unique)

  # Create vector to store RGB for each point
  cols <- id_to_color[as.character(p2)]

  # Darken foliage points
  cols = t(col2rgb(cols))
  foliage = p1 >= 1
  cols[foliage,] =   cols[foliage,] * 0.7

  # Assign to LAS
  las@data$R = cols[,1]
  las@data$G = cols[,2]
  las@data$B = cols[,3]

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





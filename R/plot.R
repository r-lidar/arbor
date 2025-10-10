#' @export
plot_anisotropy = function(las, dtm = NULL)
{
  x <- lidR::plot(las, color = "anisotropy", legend = T, breaks = "quantile")
  if (!is.null(dtm)) lidR::add_dtm3d(x, dtm)
  return(invisible(x))
}

#' @export
plot_foliage = function(las, dtm = NULL)
{
  x <- lidR::plot(las, color = "foliage", pal = foliage.colors, size = 2)
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





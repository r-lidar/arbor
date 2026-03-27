#' Convert a Quantitative Structure Forest (QSF) to an sf treemap
#'
#' Extracts global statistics from each Quantitative Structure Model (QSM)
#' contained in a Quantitative Structure Forest (QSF) object and returns
#' them as an \pkg{sf} point object. Each point corresponds to the root
#' position of an individual tree.
#'
#' @param qsf A Quantitative Structure Forest (QSF)
#'
#' @return An \pkg{sf} object with one feature per QSM. Attributes correspond
#'   to global statistics returned by \link{qsm_stats}, and geometry is
#'   defined by the root coordinates (\code{X_root}, \code{Y_root}, \code{Z_root}).
#'
#' @export
qsf_treemap = function(qsf)
{
  ans = Filter(function(x) inherits(x, "qsm"), qsf)
  ans = lapply(qsf, function(x) qsm_stats(x)$stats_global)
  ans = data.table::rbindlist(ans)
  ans = ans[,-c("X_bh", "Y_bh", "Z_bh")]
  ans = sf::st_as_sf(ans, coords = c("X_root", "Y_root", "Z_root"))
  ans
}

#' Print a Quantitative Structure Forest (QSF)
#'
#' Displays a short summary of a Quantitative Structure Forest followed by
#' the standard \code{data.table} print output.
#'
#' @param x A \code{qsm} object.
#' @param ... Unused (for S3 compatibility).
#' @method print qsf
#' @return Invisibly returns \code{x}.
#' @export
print.qsf <- function(x, ...)
{
  # number of cylinders
  n_tree <- length(x)
  n_cyl  <- sum(sapply(x, nrow))
  crs = attr(x, "crs")
  if (is.null(crs)) crs = sf::NA_crs_

  # total volume
  total_volume <- sum(sapply(x, function(x) sum(x$volume, na.rm = TRUE)))

  cat("QSF\n")
  cat(sprintf("Trees       : %d\n", n_tree))
  cat(sprintf("Cylinders   : %d\n", n_cyl))
  cat(sprintf("Volume      : %.2f m\u00b3\n", total_volume))
  cat("Coord. ref. :", crs$Name, "\n")

  invisible(x)
}

#' Print a Quantitative Structure Model (QSM)
#'
#' Displays a short summary of a Quantitative Structure Model followed by
#' the standard \code{data.table} print output.
#'
#' @param x A \code{qsm} object.
#' @param ... Unused (for S3 compatibility).
#' @method print qsm
#' @return Invisibly returns \code{x}.
#' @export
print.qsm <- function(x, ...)
{
  n_cyl <- nrow(x)
  total_volume <- 0
  height <- 0
  dbh <- 0
  if (n_cyl > 0)
  {
    total_volume <- sum(x$volume, na.rm = TRUE)
    zmin <- min(c(x$startZ, x$endZ), na.rm = TRUE)
    zmax <- max(c(x$startZ, x$endZ), na.rm = TRUE)
    height <- zmax - zmin
    dbh = qsm_dbh(x)$dbh
  }
  crs = attr(x, "crs")
  if (is.null(crs)) crs = sf::NA_crs_

  cat(sprintf("Class       : QSM\n"))
  cat(sprintf("Cylinders   : %d\n", n_cyl))
  cat(sprintf("Diameter    : %.1f cm\n", dbh*100))
  cat(sprintf("Height      : %.1f m\n", height))
  cat(sprintf("Volume      : %.2f m\u00b3\n", total_volume))
  cat("Coord. ref. :", crs$Name, "\n")

  invisible(x)
}

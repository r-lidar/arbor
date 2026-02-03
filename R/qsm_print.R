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
  cat("Quantitative Structure Model:\n")

  # number of cylinders
  n_cyl <- nrow(x)

  # total volume
  total_volume <- sum(x$volume, na.rm = TRUE)

  # tree height (robust definition)
  zmin <- min(c(x$startZ, x$endZ), na.rm = TRUE)
  zmax <- max(c(x$startZ, x$endZ), na.rm = TRUE)
  height <- zmax - zmin

  cat(sprintf("Cylinders : %d\n", n_cyl))
  cat(sprintf("Height    : %.1f m\n", height))
  cat(sprintf("Volume    : %.2f m³\n\n", total_volume))

  invisible(x)
}

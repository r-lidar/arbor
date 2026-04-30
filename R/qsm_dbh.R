#' Extract diameter at breast height (DBH) from a QSM
#'
#' Extract the diameter at breast height (DBH) from a Quantitative Structure Model
#' (QSM).
#'
#' @param qs A QSM or QSF
#' @param bh Numeric. Breast height (in meters). Default is 1.30.
#'
#' @return
#' A data.frame with:
#' \itemize{
#'   \item{dbh: DBH in meter}
#'   \item{pos: xyz position of the center}
#'   \item{normal: normal vector of the circle}
#' }
#' @examples
#' f <- system.file("extdata", "tree_qsm.laz", package="arbor")
#' tree <- lidR::readLAS(f)
#' qsm <- qsm(tree)
#' ans <- qsm_dbh(qsm)
#' \dontrun{
#' plot_semantic(tree) |> add_dbh3d(ans)
#' plot_qsm(qsm) |> add_dbh3d(ans)}
#' @export
#' @rdname dbh
qsm_dbh <- function(qs, bh = 1.30)
{
  UseMethod("qsm_dbh")
}

#' @export
qsm_dbh.qsm <- function(qs, bh = 1.30)
{
  o = qsm_dbh_cpp(qs)
  return(o)
}

#' @export
qsm_dbh.qsf <- function(qs, bh = 1.30)
{
  ans <- lapply(qs, qsm_dbh, bh = bh)
  ans <- data.table::rbindlist(ans, idcol = "treeID")
  return(ans)
}

#' Simple DBH plot function
#'
#' @param x An alignment offset returned by a `plot` function.
#' @param df the output of `qsm_dbh()`
#' @param col,lwd render parameters
#'
#' @rdname dbh
#' @export
add_dbh3d <- function(x, df, col = "red", lwd = 3)
{
  n = 32

  crossprod3 <- function(a, b)
  {
    {
      c(
        a[2] * b[3] - a[3] * b[2],
        a[3] * b[1] - a[1] * b[3],
        a[1] * b[2] - a[2] * b[1]
      )
    }
  }

  dx = x[1]
  dy = x[2]
  stopifnot(all(c("dbh", "x", "y", "z", "nx", "ny", "nz") %in% names(df)))

  for (i in seq_len(nrow(df)))
  {
    center <- c(df$x[i]-dx, df$y[i]-dy, df$z[i])
    radius <- df$dbh[i] / 2

    normal <- c(df$nx[i], df$ny[i], df$nz[i])
    normal <- normal / sqrt(sum(normal^2))

    if (abs(normal[1]) < 0.9) {
      v <- c(1, 0, 0)
    } else {
      v <- c(0, 1, 0)
    }

    # Orthonormal basis of the circle plane
    u <- crossprod3(normal, v)
    u <- u / sqrt(sum(u^2))

    w <- crossprod3(normal, u)
    w <- w / sqrt(sum(w^2))

    theta <- seq(0, 2 * pi, length.out = n)

    pts <- vapply(
      theta,
      function(t) center + radius * cos(t) * u + radius * sin(t) * w,
      numeric(3)
    )

    rgl::lines3d(x = pts[1, ], y = pts[2, ], z = pts[3, ], col = col, lwd = lwd)
  }

  return(invisible(x))
}



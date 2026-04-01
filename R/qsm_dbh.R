#' Extract diameter at breast height (DBH) from a QSM
#'
#' Extract the diameter at breast height (DBH) from a Quantitative Structure Model
#' (QSM).
#'
#' @param qs A QSM or QSF
#' @param bh Numeric. Breast height (in meters). Default is 1.30.
#'
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
#' @export
qsm_dbh <- function(qs, bh = 1.30)
{
  UseMethod("qsm_dbh")
}

#' @export
qsm_dbh.qsm <- function(qs, bh = 1.30)
{
  main_axis <- qs[qs$branch_order == 1,]
  ground_z  <- min(qs$startZ)

  len = with(main_axis, sqrt((endX - startX)^2 + (endY - startY)^2 + (endZ - startZ)^2))
  dist_to_root = cumsum(len)
  idx = which.min(abs(dist_to_root-bh))
  cyl <- main_axis[idx, ]
  if (nrow(cyl) != 1) stop("Invalid number of cylinders at breast height")

  # QSM DBH
  dbh_qsm <- round(cyl$radius * 2, 3)
  x_dbh   <- (cyl$startX + cyl$endX)/2
  y_dbh   <- (cyl$startY + cyl$endY)/2
  z_dbh   <- (cyl$startZ + cyl$endZ)/2

  # Axis geometry
  start <- c(cyl$startX, cyl$startY, cyl$startZ)
  end   <- c(cyl$endX,   cyl$endY,   cyl$endZ)

  basis <- .axis_basis(start, end)
  u <- basis$u

  o = data.frame(dbh = dbh_qsm, x = x_dbh, y = y_dbh, z = z_dbh, nx = u[1], ny = u[2], nz = u[3])
  return(o)
}

#' @export
qsm_dbh.qsf <- function(qs, bh = 1.30)
{
  ans <- lapply(qs, qsm_dbh, bh = bh)
  ans <- data.table::rbindlist(ans, idcol = "treeID")
  return(ans)
}

.axis_basis <- function(start, end)
{
  v <- end - start
  v_len <- sqrt(sum(v^2))
  if (v_len < 1e-8) stop("Degenerate cylinder axis")

  u <- v / v_len

  ref <- if (abs(u[1]) < 0.9) c(1, 0, 0) else c(0, 1, 0)

  a <- c(
    u[2]*ref[3] - u[3]*ref[2],
    u[3]*ref[1] - u[1]*ref[3],
    u[1]*ref[2] - u[2]*ref[1]
  )
  a <- a / sqrt(sum(a^2))

  b <- c(
    u[2]*a[3] - u[3]*a[2],
    u[3]*a[1] - u[1]*a[3],
    u[1]*a[2] - u[2]*a[1]
  )

  list(u = u, a = a, b = b)
}

.project_to_plane <- function(XYZ, origin, a, b)
{
  d <- sweep(XYZ, 2, origin)
  data.table::data.table(
    U = d %*% a,
    V = d %*% b,
    Z = 0
  )
}

.unproject_from_plane <- function(UV, origin, a, b)
{
  # Ensure matrix form
  UV <- as.matrix(UV)

  XYZ <- sweep(
    UV[, 1] %*% t(a) + UV[, 2] %*% t(b),
    2,
    origin,
    "+"
  )

  data.table::data.table(
    X = XYZ[, 1],
    Y = XYZ[, 2],
    Z = XYZ[, 3]
  )
}

#' @importFrom data.table :=
.extract_slice <- function(points, P_bh, u, thickness)
{
  X <- Y <- Z <- dist <- NULL
  pts <- lidR::filter_poi(points, Z >= (P_bh[3] - 0.5) &  Z <= (P_bh[3] + 0.5))

  if (nrow(pts) == 0) return(NULL)

  suppressWarnings(pts@data[, dist := (X - P_bh[1]) * u[1] + (Y - P_bh[2]) * u[2] + (Z - P_bh[3]) * u[3]])

  lidR::filter_poi(pts, abs(dist) <= thickness / 2)
}

.circle_world <- function(center_uv, radius, P_bh, a, b, n = 200)
{
  t <- seq(0, 2*pi, length.out = n)

  data.table::data.table(
    X = P_bh[1] +
      center_uv[1] * a[1] + center_uv[2] * b[1] +
      radius * cos(t) * a[1] + radius * sin(t) * b[1],

    Y = P_bh[2] +
      center_uv[1] * a[2] + center_uv[2] * b[2] +
      radius * cos(t) * a[2] + radius * sin(t) * b[2],

    Z = P_bh[3] +
      center_uv[1] * a[3] + center_uv[2] * b[3] +
      radius * cos(t) * a[3] + radius * sin(t) * b[3]
  )
}



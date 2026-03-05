#' Estimate diameter at breast height (DBH) from a QSM and point cloud slice
#'
#' Computes diameter at breast height (DBH) from a Quantitative Structure Model
#' (QSM). When the corresponding point cloud is provided, an additional DBH estimate is
#' obtained by extracting a slice of the point cloud at breast height and fitting a
#' circle using RANSAC. This is similar to the method used for computing the QSM but the
#' QSM radii are smoothed with polynomial fitting. If \code{tree} is \code{NULL}, only
#' the QSM-based DBH is returned.
#'
#' @param qsm A \code{data.frame} or \code{data.table} describing a QSM.
#' @param tree Optional \code{LAS} object containing the point cloud of the tree.
#'   If provided, DBH is also estimated from a point cloud slice.
#' @param slice_thickness Numeric. Thickness (in meters) of the slice extracted
#' around breast height.
#' @param bh Numeric. Breast height (in meters). Default is 1.37.
#' @param display Logical. If \code{TRUE}, diagnostic plots comparing QSM and
#'   slice-based DBH estimates are produced.
#'
#' @details
#' Breast height is computed relative to the minimum \code{startZ} value in the
#' QSM. If the QSM is not propagated to the ground this value may be erroneous.
#' When \code{tree} is supplied, points are projected onto the plane orthogonal
#' to the main axis at breast height, and a circle is fitted using a RANSAC-based
#' method.\cr\cr
#' Diagnostic plots include vertical and horizontal views of the slice, QSM
#' cylinder, and fitted circles.
#'
#' @return
#' A named nested list with:
#' \itemize{
#'   \item{qsm: Derived from the QSM cylinder}
#'   \item{slice: Derived from the point cloud slice.}
#'   \item{average: Mean of QSM and slice-based.}
#' }
#' Each item contains:
#' \itemize{
#'   \item{dbh: DBH in meter}
#'   \item{pos: xyz position of the center}
#'   \item{normal: normal vector of the circle}
#' }
#' @examples
#' f <- system.file("extdata", "tree_qsm.laz", package="arbor")
#' tree <- lidR::readLAS(f)
#' qsm <- qsm(tree)
#' ans <- qsm_dbh(qsm, tree, display = TRUE)
#' @export
qsm_dbh <- function(qsm, tree = NULL, slice_thickness = 0.1, bh = 1.30, display = FALSE)
{
  X <- Y <- Z <- branch_order <- startZ <- endZ <- . <- NULL
  data.table::setDT(qsm)

  main_axis <- qsm[branch_order == 1]
  ground_z  <- min(qsm$startZ)

  len = with(qsm, sqrt((endX - startX)^2 + (endY - startY)^2 +(endZ - startZ)^2))
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
  u <- basis$u; a <- basis$a; b <- basis$b

  oqsm = list(dbh = dbh_qsm, pos = c(x = x_dbh, y = y_dbh, z = z_dbh), normal = u)
  oslice = list(dbh = NA_real_, pos = c(x = NA_real_, y = NA_real_, z = NA_real_), normal =  u)
  oavg = list(dbh = dbh_qsm,  pos = c(x = x_dbh, y = y_dbh, z = z_dbh), normal = u)

  if (is.null(tree)) return(list(qsm = oqsm, slice = oslice, average = oavg))

  if (!methods::is(tree, "LAS")) stop("'tree' must be a LAS object")

  tree  <- filter_tree(tree)

  # Interpolate BH origin
  t <- (bh+ground_z - start[3]) / (end[3] - start[3])
  P_bh <- start + t * (end - start)

  # Slice extraction
  slice <- .extract_slice(tree, P_bh, u, slice_thickness)

  if (is.null(slice))
    return(list(qsm = oqsm, slice = oslice, average = oavg))

  # This should be computed reprojected!!
  d_center = sqrt((slice$X-x_dbh)^2 + (slice$Y-y_dbh)^2)
  keep = d_center < dbh_qsm
  slice = slice[keep]

  XYZ <- as.matrix(slice@data[, .(X, Y, Z)])
  UV  <- .project_to_plane(XYZ, P_bh, a, b)

  if (nrow(UV) < 5)
  {
    warning("Impossible to fit a circle on a slice with less that 5 points. Returned QSM's DBH only.")
    return(list(qsm = oqsm, slice = oslice, average = oavg))
  }

  # Circle fit
  fit <- ransac_circle(as.matrix(UV), 1000, 0.02)

  if (fit$covered_arc_degree < 90 | fit$percentage_inlier * 100 < 0.5)
  {
    warning("Impossible to fit a valid circle with RANSAC. Returned QSM's DBH only.")
    return(list(qsm = oqsm, slice = oslice, average = oavg))
  }

  r = fit$radius

  ratio = (abs((2*r) - dbh_qsm))/dbh_qsm
  if (ratio > 0.5)
  {
    warning("Too different radii between QSM slice measurement. Returned QSM's DBH only.")
    return(list(qsm = oqsm, slice = oslice, average = oavg))
  }

  dbh_slice <- unname(round(r * 2, 3))
  px_dbh2 = unname(fit$center_x)
  py_dbh2 = unname(fit$center_y)
  xy_dbh2 = matrix(c(px_dbh2, py_dbh2, 0), 1, 3)
  xy_dbh2 = .unproject_from_plane(xy_dbh2, P_bh, a, b)
  x_dbh2 = as.numeric(xy_dbh2[1,1])
  y_dbh2 = as.numeric(xy_dbh2[1,2])
  z_dbh2 = as.numeric(xy_dbh2[1,3])

  oslice = list(dbh = dbh_slice,
                pos = c(x = x_dbh2, y = y_dbh2, z = z_dbh2),
                normal = u)
  oavg = list(dbh = (dbh_slice+dbh_qsm)/2,
              pos = c(x = (x_dbh2+x_dbh)/2, y = (y_dbh2+y_dbh)/2, z = (z_dbh2+z_dbh)/2),
              normal = u)

  # ---- PLOTS ----
  if (display)
  {
    opar <- graphics::par(
      mfrow = c(3, 2),
      mar   = c(4, 4, 2.5, 1),
      mgp   = c(2.2, 0.7, 0),
      tcl   = -0.3,
      cex.lab  = 0.95,
      cex.axis = 0.85
    )
    on.exit(graphics::par(opar))

    butt = lidR::filter_poi(tree, Z < bh+ground_z + 1)

    # World diagnostic
    circle_xyz1 <- .circle_world(c(0, 0), dbh_qsm/2, P_bh, a, b)
    circle_xyz2 <- .circle_world(c(px_dbh2, py_dbh2), dbh_slice/2, P_bh, a, b)

    graphics::plot(butt@data$X, butt@data$Z, asp = 1, pch = 19, cex = 0.3,  xlab = "X", ylab = "Z", main = "", ylim = c(ground_z-0.1, max(butt$Z)+0.1))
    graphics::points(slice@data$X, slice@data$Z, asp = 1, pch = 19, cex = 0.4, col = "purple")
    graphics::lines(circle_xyz1$X, circle_xyz1$Z, col = "green", lwd = 2)
    graphics::lines(circle_xyz2$X, circle_xyz2$Z, col = "red", lwd = 2)
    graphics::abline(h = ground_z, col = "black", lty = 1, lwd = 2)

    graphics::plot(butt@data$Y, butt@data$Z, asp = 1, pch = 19, cex = 0.3,  xlab = "Y", ylab = "Z", main = "", ylim = c(ground_z-0.1, max(butt$Z)+0.1))
    graphics::points(slice@data$Y, slice@data$Z, asp = 1, pch = 19, cex = 0.4, col = "purple")
    graphics::lines(circle_xyz1$Y, circle_xyz1$Z, col = "green", lwd = 2)
    graphics::lines(circle_xyz2$Y, circle_xyz2$Z, col = "red", lwd = 2)
    graphics::abline(h = ground_z, col = "black", lty = 1, lwd = 2)

    butt = lidR::filter_poi(tree, Z > bh+ground_z - 0.2, Z < bh+ground_z + 0.2)

    graphics::plot(butt@data$X, butt@data$Z, asp = 1, pch = 19, cex = 0.2,  xlab = "X", ylab = "Z", main = "")
    graphics::points(slice@data$X, slice@data$Z, asp = 1, pch = 19, cex = 0.5, col = "purple")
    graphics::lines(circle_xyz1$X, circle_xyz1$Z, col = "green", lwd = 2)
    graphics::lines(circle_xyz2$X, circle_xyz2$Z, col = "red", lwd = 2)

    graphics::plot(butt@data$Y, butt@data$Z, asp = 1, pch = 19, cex = 0.2,  xlab = "Y", ylab = "Z", main = "")
    graphics::points(slice@data$Y, slice@data$Z, asp = 1, pch = 19, cex = 0.5, col = "purple")
    graphics::lines(circle_xyz1$Y, circle_xyz1$Z, col = "green", lwd = 2)
    graphics::lines(circle_xyz2$Y, circle_xyz2$Z, col = "red", lwd = 2)
    graphics::abline(h = ground_z, col = "black", lty = 1, lwd = 2)


    graphics::plot(slice@data$X, slice@data$Y, asp = 1, pch = 19, cex = 0.4, col = "purple", xlab = "X", ylab = "Y")
    graphics::lines(circle_xyz1$X, circle_xyz1$Y, col = "green", lwd = 2)
    graphics::lines(circle_xyz2$X, circle_xyz2$Y, col = "red", lwd = 2)

    graphics::plot.new()

    graphics::legend(
      "center",
      legend = c(
        "Slice points",
        "Extracted from QSM",
        "Computed from slice"
      ),
      col = c("purple", "green", "red"),
      pch = c(19, NA, NA),
      lty = c(NA, 1, 1),
      lwd = c(NA, 2, 2),
      pt.cex = 1.2,
      cex = 1.1,
      bty = "n"
    )
  }

  return(list(qsm = oqsm, slice = oslice, average = oavg))
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



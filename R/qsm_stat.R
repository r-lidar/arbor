#' Compute summary statistics from a QSM
#'
#' Computes structural statistics from a Quantitative Structure Model (QSM),
#' including volume distributions by radius and branch order, stem taper,
#' and global tree metrics such as DBH, height, and total volume. Optionally,
#' statistics can be complemented with information derived from a LAS point cloud
#' of the tree and visual diagnostics can be displayed.
#'
#' @param qsm A QSM or QSF produced by \link{qsm} or \link{qsf}.
#' @param display Logical. If `TRUE`, diagnostic plots are displayed, including
#'   volume distributions, stem profile, and vertical projections.
#' @param ... passed to \link{qsm_dbh}.
#'
#' @details
#' The function computes:
#' \itemize{
#'   \item Volume and cylinder counts by radius class (2 cm bins).
#'   \item Volume distribution and cumulative percentage by branch order.
#'   \item Stem taper profile from the main axis.
#'   \item Global tree statistics such as DBH, total volume, and height.
#' }
#' When a `LAS` object is supplied, metrics are derived directly
#' from the point cloud and QSM rather than the QSM alone.
#'
#' @return A named list with the following elements:
#' \itemize{
#'   \item{stats_by_order: A `data.table` with volume and cumulative percentage by branch order.}
#'   \item{stats_by_radius: A `data.table` with total volume and cylinder count by radius class.}
#'   \item{stats_global: A `data.frame` containing global tree metrics (coordinates, DBH, volume, height).}
#'   \item{stem_taper: A `data.frame` describing stem diameter as a function of distance to root.}
#' }
#'
#' @seealso
#' \link{qsm_dbh}, \link{qsm}
#' @examples
#' f <- system.file("extdata", "tree_qsm.laz", package="arbor")
#' tree <- lidR::readLAS(f)
#' qsm <- qsm(tree)
#' ans <- qsm_stats(qsm, tree, display = TRUE)
#' @export
qsm_stats <- function(qs, ..., display = FALSE)
{
  UseMethod("qsm_stats")
}

#' @export
qsm_stats.qsm <- function(qs, ..., display = FALSE)
{
  if (nrow(qs) == 0) return(NULL)

  if (!is.data.frame(qs)) stop("'qs' must be a data.frame or data.table")

  required_cols <- c("radius", "branch_order", "axis_ID", "startX", "startY", "startZ", "endZ")
  missing_cols <- setdiff(required_cols, names(qs))
  if (length(missing_cols) > 0) stop("Missing required QSM columns: ", paste(missing_cols, collapse = ", "))

  if (!is.logical(display) || length(display) != 1)  stop("'display' must be a single logical value (TRUE/FALSE)")

  .N <- radius <- radius_bin <- branch_order <- axis_ID <- dist_to_root <-. <- NULL

  dt = data.table::copy(qs)
  dt$radius_bin = cut(dt$radius*100, breaks = seq(0, max(dt$radius*100)+2, 2))

  stats_by_radius <- dt[, .(volume = round(sum(volume),3), N = .N), by = .(radius = radius_bin)]
  data.table::setDT(stats_by_radius)
  data.table::setorder(stats_by_radius, radius)

  stats_by_order <- dt[, .(volume = round(sum(volume), 3)), by = branch_order]
  data.table::setDT(stats_by_order)
  data.table::setorder(stats_by_order, branch_order)
  stats_by_order$percentage = round(cumsum(stats_by_order$volume)/sum(stats_by_order$volume)*100,1)

  main_axis = dt[axis_ID == 1]
  sl = main_axis$subtree_length
  dist_to_root = max(sl)-sl
  r = main_axis$radius
  stem_profile = data.frame(dist_to_root = dist_to_root, diametre = r*2)
  data.table::setDT(stem_profile)
  data.table::setorder(stem_profile, dist_to_root)

  dbh = qsm_dbh(qs, ...)
  X_bh = dbh$x
  Y_bh = dbh$y
  Z_bh = dbh$z

  root = qs[axis_ID == 1][1]
  X_root = root$startX
  Y_root = root$startY
  Z_root = root$startZ

  qs = qsm_volume(qs)
  V   = round(sum(qs$volume), 3)
  H   = round(max(qs$endZ) - Z_root, 3)
  DBH = dbh$dbh

  stats_global = data.frame(X_root = X_root,
                            Y_root = Y_root,
                            Z_root = Z_root,
                            X_bh = X_bh,
                            Y_bh = Y_bh,
                            Z_bh = Z_bh,
                            DBH = DBH,
                            V = V,
                            H = H)

  if (display) .plot_stats(stats_by_radius, stats_by_order, stats_global, stem_profile, NULL)

  return(list(stats_by_order = stats_by_order,
              stats_by_radius = stats_by_radius,
              stats_global = stats_global,
              stem_taper = stem_profile))

}

#' @export
qsm_stats.qsf <- function(qs, ..., display = FALSE)
{
  if (!is.list(qs)) stop("'qsm' must be a list of qsm objects")

  res <- lapply(qs, function(x)
  {
    out <- qsm_stats(x, ..., display = FALSE)
    if (is.null(out)) return(NULL)
    out$stats_global
  })

  res <- Filter(Negate(is.null), res)
  if (length(res) == 0) return(NULL)

  res <- data.table::rbindlist(res, idcol = "treeID")

  if (display) {
    message("Display not supported for qsf (global stats only).")
  }

  return(res)
}

.plot_stats = function(stats_by_radius, stats_by_order, stats_global, stem_profile, tree)
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

  # 1) Volume by radius
  graphics::barplot(
    stats_by_radius$volume,
    names.arg = stats_by_radius$radius,
    col       = "grey70",
    border    = "grey30",
    xlab      = "Radius (cm)",
    ylab      = expression("Total volume (m"^3*")"),
    main      = "Volume by radius"
  )
  graphics::box()

  # 2) Number of cylinders by radius
  graphics::barplot(
    stats_by_radius$N,
    names.arg = stats_by_radius$radius,
    col       = "grey80",
    border    = "grey30",
    xlab      = "Radius (cm)",
    ylab      = "Number of cylinders",
    main      = "Cylinder count by radius"
  )
  graphics::box()

  # 3) Volume by branch order
  graphics::barplot(
    stats_by_order$volume,
    names.arg = stats_by_order$branch_order,
    col       = "grey70",
    border    = "grey30",
    xlab      = "Branch order",
    ylab      = expression("Total volume (m"^3*")"),
    main      = "Volume by branch order"
  )
  graphics::box()

  # 4) Volume percentage by branch order
  graphics::barplot(
    stats_by_order$percentage,
    names.arg = stats_by_order$branch_order,
    col       = "grey85",
    border    = "grey30",
    xlab      = "Branch order",
    ylab      = "Total volume (%)",
    main      = "Relative volume contribution"
  )
  graphics::box()

  # 5) Stem profile
  graphics::plot(
    stem_profile$diametre*100, stem_profile$dist_to_root,
    pch  = 19,
    cex  = 0.45,
    col  = grDevices::adjustcolor("black", 0.5),
    xlab = "Diameter (cm)",
    ylab = "Distance to root (m)",
    main = "Stem profile"
  )
  graphics::grid(col = "grey85", lty = 1)

  # 6) Tree vertical structure
  if (!is.null(tree))
  {
    graphics::plot(
      tree$X, tree$Z - stats_global$Z_root,
      pch  = 19,
      cex  = 0.05,
      col  = grDevices::adjustcolor("black", 0.4),
      asp  = 1,
      ylim = c(- 0.5, stats_global$H + 0.5),
      xlab = "Horizontal distance (m)",
      ylab = "Height (m)",
      main = "Tree vertical structure"
    )

    graphics::abline(h = 0, lwd = 2)
    graphics::abline(h = stats_global$H, lty = 3)
    graphics::abline(h = stats_global$Z_bh-stats_global$Z_root, lty = 3, col = "steelblue")
    graphics::box()
  }
}

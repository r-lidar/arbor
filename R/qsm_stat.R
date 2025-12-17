#' Compute summary statistics from a QSM
#'
#' Computes structural statistics from a Quantitative Structure Model (QSM),
#' including volume distributions by radius and branch order, stem taper,
#' and global tree metrics such as DBH, height, and total volume. Optionally,
#' statistics can be complemented with information derived from a LAS point cloud
#' of the tree and visual diagnostics can be displayed.
#'
#' @param qsm A `data.table` or `data.frame` representing a QSM, typically
#'   produced by \link{qsm}.
#' @param tree Optional. A `LAS` object containing the tree point cloud that generated the QSM
#'   If provided, additional statistics are computed.
#' @param display Logical. If `TRUE`, diagnostic plots are displayed, including
#'   volume distributions, stem profile, and vertical projections.
#'
#' @details
#' The function computes:
#' \itemize{
#'   \item Volume and cylinder counts by radius class (2 cm bins).
#'   \item Volume distribution and cumulative percentage by branch order.
#'   \item Stem taper profile from the main axis.
#'   \item Global tree statistics such as DBH, total volume, and height.
#' }\cr\cr
#' When a `LAS` object is supplied, metrics are derived directly
#' from the point cloud and QSM rather than the QSM alone.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{stats_by_order}{A `data.table` with volume and cumulative percentage by branch order.}
#'   \item{stats_by_radius}{A `data.table` with total volume and cylinder count by radius class.}
#'   \item{stats_global}{A `data.frame` containing global tree metrics (coordinates, DBH, volume, height).}
#'   \item{stem_taper}{A `data.frame` describing stem diameter as a function of distance to root.}
#' }
#'
#' @seealso
#' \code{\link{qsm_dbh}}, \code{\link{qsm_qsm}}
#'
#' @export
qsm_stats <- function(qsm, tree = NULL, display = FALSE)
{
  if (missing(qsm)) stop("'qsm' is missing")
  if (!is.data.frame(qsm)) stop("'qsm' must be a data.frame or data.table")

  required_cols <- c("radius", "branch_order", "axis_ID", "startX", "startY", "startZ", "endZ")
  missing_cols <- setdiff(required_cols, names(qsm))
  if (length(missing_cols) > 0) stop("Missing required QSM columns: ", paste(missing_cols, collapse = ", "))

  if (!is.null(tree)) {
    if (!methods::is(tree, "LAS"))
      stop("'tree' must be a LAS object")
  }

  if (!is.logical(display) || length(display) != 1)  stop("'display' must be a single logical value (TRUE/FALSE)")

  .N <- radius <- radius_bin <- branch_order <- axis_ID <- dist_to_root <-. <- NULL

  dt = data.table::copy(qsm)
  dt$radius_bin = cut(dt$radius*100, breaks = seq(0, max(dt$radius*100)+2, 2))

  stats_by_radius <- dt[, .(volume = round(sum(volume),3), N = .N), by = .(radius = radius_bin)]
  data.table::setorder(stats_by_radius, radius)

  stats_by_order <- dt[, .(volume = round(sum(volume), 3)), by = branch_order]
  data.table::setorder(stats_by_order, branch_order)
  stats_by_order$percentage = round(cumsum(stats_by_order$volume)/sum(stats_by_order$volume)*100,1)

  main_axis = dt[axis_ID == 1]
  sl = main_axis$subtree_length
  dist_to_root = max(sl)-sl
  r = main_axis$radius
  stem_profile = data.frame(dist_to_root = dist_to_root, diametre = r*2)
  data.table::setorder(stem_profile, dist_to_root)

  dbh = qsm_dbh(qsm, tree, display = display)
  X_bh = dbh$average$x
  Y_bh = dbh$average$y
  Z_bh = dbh$average$z

  root = qsm[axis_ID == 1][1]
  X_root = root$startX
  Y_root = root$startY
  Z_root = root$startZ

  qsm = qsm_volume(qsm)
  V   = round(sum(qsm$volume), 3)
  H   = round(max(qsm$endZ) - Z_root, 3)

  stats_global = data.frame(X_root = X_root,
                            Y_root = Y_root,
                            Z_root = Z_root,
                            X_bh = Y_bh,
                            Y_bh = Y_bh,
                            Z_bh = Z_bh,
                            DBH = dbh$average$dbh,
                            V = V,
                            H = H)

  if (!is.null(tree))
  {
    H = max(tree$Z)
    q98 = quantile(tree$Z, probs = 0.98)
    stats_global$H = H - Z_root
    stats_global$q98 = q98- Z_root
  }

  if (display)
  {
    opar = graphics::par(mfrow = c(3, 2))
    on.exit(graphics::par(opar))

    graphics::barplot(stats_by_radius$volume, names.arg = stats_by_radius$radius, xlab = "Radius (cm)", ylab = "Total volume (m\u00B3)")
    graphics::barplot(stats_by_radius$N, names.arg = stats_by_radius$radius, xlab = "Radius(cm)", ylab = "Num. cylinders")

    graphics::barplot(stats_by_order$volume, names.arg = stats_by_order$branch_order, xlab = "Branch order", ylab = "Total volume (m\u00B3)")
    graphics::barplot(stats_by_order$percentage, names.arg = stats_by_order$branch_order, xlab = "Branch order", ylab = "Total volume (%)")

    graphics::plot(r*200, dist_to_root, xlab = "Diameter (cm)", ylab = "Distance to root (m)", main = "Stem profile", pch = 19, cex = 0.5)

    graphics::plot(tree$X, tree$Z, pch = 19, cex = 0.05, asp = 1, ylim = c(Z_root-0.5, stats_global$H+Z_root+0.5))
    graphics::abline(h = Z_root, lty = 1, lwd = 2)
    graphics::abline(h = stats_global$H+Z_root, lty = 3)
    graphics::abline(h = stats_global$Z_bh, lty = 3, col = "blue")
  }

  return(list(stats_by_order = stats_by_order,
              stats_by_radius = stats_by_radius,
              stats_global = stats_global,
              stem_taper = stem_profile))

}

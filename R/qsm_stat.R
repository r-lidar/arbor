qsm_stat = function(qsm)
{
  .N <- radius_bin <- branch_order <- axis_ID <- . <- NULL

  opar = graphics::par(mfrow = c(3, 2))
  on.exit(graphics::par(opar))

  dt = data.table::copy(qsm)
  dt$radius_bin = cut(dt$radius*100, breaks = seq(0, max(dt$radius*100)+2, 2))

  out <- dt[, .(total_volume = sum(volume), N = .N), by = radius_bin]
  data.table::setorder(out, radius_bin)
  graphics::barplot(out$total_volume, names.arg = out$radius_bin, xlab = "Radius (cm)", ylab = "Total volume (m\u00B3)")
  graphics::barplot(out$N, names.arg = out$radius_bin, xlab = "Radius(cm)", ylab = "Num. cylinders")

  out <- dt[, .(total_volume = sum(volume)), by = branch_order]
  data.table::setorder(out, branch_order)
  graphics::barplot(out$total_volume, names.arg = out$branch_order, xlab = "Branch order", ylab = "Total volume (m\u00B3)")

  out$cum_volume = cumsum(out$total_volume)/sum(out$total_volume)*100
  graphics::barplot(out$cum_volume, names.arg = out$branch_order, xlab = "Branch order", ylab = "Total volume (%)")

  main_axis = dt[axis_ID == 1]
  sl = main_axis$subtree_length
  dist_to_root = max(sl)-sl
  r = main_axis$radius*100
  plot(r, dist_to_root, xlab = "Radius (cm)", ylab = "Distance to root (m)", main = "Stem profile", pch = 19, cex = 0.5)

}

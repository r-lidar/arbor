qsm_prolongation <- function(qsm, d, L = 0.1)
{
  cum_length <- cyl_length <- NULL

  if (d <= 0) return(qsm)

  qsm$cyl_length = sqrt((qsm$endX-qsm$startX)^2+(qsm$endY-qsm$startY)^2+(qsm$endZ-qsm$startZ)^2)

  main_axis <- qsm[qsm$axis_ID == 1]
  root <- main_axis[main_axis$parent_ID == 0, ]
  main_axis[, cum_length := cumsum(cyl_length)]
  dt_top <- main_axis[cum_length <= 0.1*max(main_axis$cum_length)]

  if (sum(dt_top$cyl_length) < 0.3)
  {
    next_cyl <- main_axis[nrow(dt_top) + 1]
    dt_top <- rbind(dt_top, next_cyl)
  }

  dx <- dt_top$endX[nrow(dt_top)] - dt_top$startX[1]
  dy <- dt_top$endY[nrow(dt_top)] - dt_top$startY[1]
  dz <- dt_top$endZ[nrow(dt_top)] - dt_top$startZ[1]
  norm <- sqrt(dx^2 + dy^2 + dz^2)
  orientation <- c(dx, dy, dz) / norm
  d <- d / orientation[3]

  end_point <- c(root$startX, root$startY, root$startZ)
  start_point <- end_point - orientation * d

  # Number of segments
  n_segments <- max(1, ceiling(d / L))
  actual_L <- d / n_segments

  # Generate subdivided cylinders
  new_cyls <- vector("list", n_segments)
  for (i in seq_len(n_segments))
  {
    frac1 <- (i - 1) / n_segments
    frac2 <- i / n_segments

    pt1 <- end_point - orientation * d * frac1
    pt2 <- end_point - orientation * d * frac2

    cyl_id <- -i
    parent_id <- if (i == 1) root$cyl_ID else - (i - 1)

    new_cyls[[i]] <- data.table::data.table(
      startX = pt1[1],
      startY = pt1[2],
      startZ = pt1[3],
      endX = pt2[1],
      endY = pt2[2],
      endZ = pt2[3],
      cyl_ID = cyl_id,
      parent_ID = parent_id,
      axis_ID = 1,
      subtree_length = root$subtree_length + d - actual_L * (n_segments - i + 1),
      branch_order = 1,
      cyl_length = sqrt(sum((pt2 - pt1)^2))
    )
  }

  qsm <- data.table::rbindlist(c(new_cyls, list(qsm)), use.names = TRUE, fill = TRUE)
  qsm$cyl_length = NULL
  return(qsm)
}

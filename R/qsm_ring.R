qsm_ring = function(qsm)
{
  target_dist <- 1.3
  trunk = qsm[axis_ID == 1]
  len = with(trunk, sqrt((endX - startX)^2 + (endY - startY)^2 +(endZ - startZ)^2))
  dist_to_root = cumsum(len)
  idx = which.min(abs(dist_to_root-target_dist))
  target_row <- trunk[idx, ]

  current_len <- with(target_row, sqrt((endX - startX)^2 + (endY - startY)^2 + (endZ - startZ)^2))

  unit_v_x <- (target_row$endX - target_row$startX) / current_len
  unit_v_y <- (target_row$endY - target_row$startY) / current_len
  unit_v_z <- (target_row$endZ - target_row$startZ) / current_len

  mid_x <- (target_row$startX + target_row$endX) / 2
  mid_y <- (target_row$startY + target_row$endY) / 2
  mid_z <- (target_row$startZ + target_row$endZ) / 2

  new_len <- 0.03
  half_len <- new_len / 2

  target_row$startX <- mid_x - (half_len * unit_v_x)
  target_row$startY <- mid_y - (half_len * unit_v_y)
  target_row$startZ <- mid_z - (half_len * unit_v_z)

  target_row$endX   <- mid_x + (half_len * unit_v_x)
  target_row$endY   <- mid_y + (half_len * unit_v_y)
  target_row$endZ   <- mid_z + (half_len * unit_v_z)
  target_row$radius <- target_row$radius + 0.01
  target_row
}

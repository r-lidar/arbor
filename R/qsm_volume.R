qsm_volume = function(qsm)
{
  l = sqrt((qsm$endX-qsm$startX)^2+(qsm$endY-qsm$startY)^2+(qsm$endZ-qsm$startZ)^2)
  qsm$volume = pi*qsm$radius^2*l
  return(qsm)
}

qsm_stem = function(qsm, min_radius = 0.045)
{
  qsm[qsm$branch_order == 1 & qsm$radius >= min_radius]
}

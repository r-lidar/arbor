qsm_simplify <- function(qsm, max_length = 0.3)
{
  qsm$cyl_length = sqrt((qsm$startX - qsm$endX)^2 + (qsm$startY- qsm$endY)^2 + (qsm$startZ - qsm$endZ)^2)
  qsm = qsm_simplify_cpp(qsm, max_length)
  qsm = qsm_topology(qsm)
  qsm = qsm_architecture(qsm)
  qsm
}

qsm_topology = function(skeleton)
{
  logger("Computing qsm topology")

  skeleton$cyl_ID = 1:nrow(skeleton)
  skeleton$parent_ID = 0
  skeleton = qsm_topology_cpp(skeleton)

  n_root = sum(skeleton$parent_ID == 0)

  # No root is a bug
  if (n_root == 0)
    stop("Internal error: no root found")

  # Two roots are rare but possible if forking at root
  # Add a 1 mm cylinder to force one root
  if (n_root > 1)
  {
    xyz = skeleton[which(skeleton$parent_ID == 0)[1],]
    xyz$endX = xyz$startX
    xyz$endY = xyz$startY
    xyz$endZ = xyz$startZ
    xyz$startZ = xyz$startZ - 0.001
    skeleton = rbind(xyz, skeleton)
    skeleton$cyl_ID = 1:nrow(skeleton)
    skeleton$parent_ID = 0
    skeleton = qsm_topology_cpp(skeleton)
  }

  data.table::setDT(skeleton)

  return(skeleton)
}

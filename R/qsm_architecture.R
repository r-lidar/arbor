qsm_architecture <- function(qsm, node_id = 0)
{
  tmp = qsm_length(qsm)
  ans = cpp_compute_architecture(tmp)
  qsm$axis_ID = ans$axis_ID
  qsm$cyl_length = tmp$length
  qsm$subtree_length = ans$subtree_length
  qsm$branch_order = ans$branching_order
  return(qsm)
}

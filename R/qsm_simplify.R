qsm_simplify <- function(qsm, max_length = 0.3)
{
  qsm = qsm_simplify_cpp(qsm, max_length)
  qsm = qsm_topology(qsm)
  qsm = qsm_architecture(qsm)
  qsm
}

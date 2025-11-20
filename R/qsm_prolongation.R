qsm_prolongation <- function(qsm, d, L = 0.1)
{
  qsm = qsm_prolongation_cpp(qsm, d, L)
  data.table::setDT(qsm)
  return(qsm[])
}

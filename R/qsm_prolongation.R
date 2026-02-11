qsm_prolongation <- function(qsm, d, L = 0.1)
{
  logger("Prolongation")
  qsm = qsm_prolongation_cpp(qsm, d, L)
  data.table::setDT(qsm)
  return(qsm[])
}

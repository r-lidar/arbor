qsm_architecture <- function(qsm)
{
  logger("Building architecture")

  qsm <- qsm_architecture_cpp(qsm)
  data.table::setDT(qsm)

  return(qsm[])
}

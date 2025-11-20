qsm_architecture <- function(qsm)
{
  cat("Building architecture...") ; ti = tic()

  qsm <- qsm_architecture_cpp(qsm)
  data.table::setDT(qsm)

  toc(ti)

  return(qsm[])
}

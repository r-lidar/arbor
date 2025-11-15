qsm_architecture <- function(qsm)
{
  cat("Building architecture...") ; ti = tic()

  qsm <- qsm_architecture_cpp(qsm)
  qsm <- qsm_smooth(qsm, niter = 1)
  qsm <- qsm_architecture_cpp(qsm)
  data.table::setDT(qsm)

  toc(ti)

  return(qsm)
}

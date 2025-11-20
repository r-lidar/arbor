qsm_smooth = function(qsm, niter = 2)
{
  ans = qsm_smooth_cpp(qsm, niter = niter)
  data.table::setDT(ans)
  return(ans)
}

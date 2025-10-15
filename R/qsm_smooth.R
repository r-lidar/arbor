qsm_smooth = function(qsm, niter = 2)
{
  ans = cpp_smooth_skeleton(qsm, niter = niter)
  ans = data.table::as.data.table(ans)
  qsm[,1:6] = ans
  return(qsm)
}

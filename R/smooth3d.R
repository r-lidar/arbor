smooth3d = function(las, radius = 0.05, weight = NULL, ncpu = 8, progress = TRUE)
{
  ti = Sys.time()

  if (is.function(radius)) radius = radius(las$Z)
  if (length(radius) == 1L) radius = rep(radius, lidR::npoints(las))
  stopifnot(length(radius) == lidR::npoints(las))
  stopifnot(all(radius > 0))
  if (is.null(weight)) weight = rep(1, lidR::npoints(las))
  stopifnot(length(weight) == lidR::npoints(las), radius > 0)

  z = cpp_smooth3d(las, radius = radius, weight = weight, ncpu = ncpu, progress, FALSE)

  tf = Sys.time()
  print(difftime(tf,ti))

  las$X = z$X
  las$Y = z$Y
  las$Z = z$Z
  las
}

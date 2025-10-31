make_master_seed = function(las)
{
  x <- mean(las$X)
  y <- mean(las$Y)
  z <- min(las$Z)-1
  seed <- data.frame(X= x, Y = y, Z = z, anisotropy = 1, pointID = 0)
  lidR::quantize(seed[["X"]], 0.01, las@header[["X offset"]])
  lidR::quantize(seed[["Y"]], 0.01, las@header[["Y offset"]])
  lidR::quantize(seed[["Z"]], 0.01, las@header[["Z offset"]])
  header <- rlas::header_create(seed)
  seed   <- suppressWarnings(lidR::LAS(seed, header))
  seed
}

evaluate_penalty = function(params)
{
  penalty = params$path_finder$angle_penalty(0:180)
  if (length(penalty) != 181) stop("Invalid penalty function")
  params$path_finder$penalty = penalty
  params
}


compute_network = function(data, query, k = 5)
{
  if (missing(query))
  {
    nn = lidR::knn(data, k = k)
    n = lidR::npoints(data)
  }
  else
  {
    nn = lidR::knnx(data, query, k = k)
    n = lidR::npoints(query)
  }

  from <- rep(1:n, each = k)
  to <- as.vector(t(nn$nn.index))
  cost <- as.vector(t(nn$nn.dist))
  edges <- data.frame(from, to, cost)
  edges
}

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



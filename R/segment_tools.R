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

decimate_translate <- function(las, params = default_parameters)
{
  decimated <- NULL

  if (!"decimated" %in% names(las))
    dec <- barycentric_predecimation(las, params)
  else
    dec <- lidR::filter_poi(las, decimated == TRUE)

  x_translation <- round(mean(las$X))
  y_translation <- round(mean(las$Y))
  dec@data$X <- dec@data$X - x_translation
  dec@data$Y <- dec@data$Y - y_translation
  dec@header@VLR$Extra_Bytes <- NULL
  return(list(dec = dec, x_translation = x_translation, y_translation = y_translation))
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

compute_point_network <- function(dec, k, max_gap = 1, wood_mask = NULL, cost_factors = NULL, power = 3, downward = FALSE)
{
  if (!is.null(wood_mask))
  {
    stopifnot(
      !is.null(cost_factors))
  }

  if (!is.null(cost_factors))
  {
    stopifnot(
      !is.null(cost_factors$wood2wood),
      !is.null(cost_factors$wood2leaf),
      !is.null(cost_factors$leaf2leaf))
  }

  net <- compute_network(dec, k = k)

  # Gaps have an infinite cost
  net$cost[net$cost > max_gap] <- Inf

  # It is more expensive to move in large steps
  net$cost <- net$cost^power

  # It is more expensive to move downward. The path finder starts from
  # a seed below the ground. If should reach any target points by moving upward
  # preferentially. If it moves downward it is not following a tree. Maybe a
  # branch bending. It is ok but expensive.
  X1 <- dec$X[net$from]; X2 <- dec$X[net$to]
  Y1 <- dec$Y[net$from]; Y2 <- dec$Y[net$to]
  Z1 <- dec$Z[net$from]; Z2 <- dec$Z[net$to]
  dz <- Z1-Z2
  magnitude <- sqrt((X1-X2)^2 + (Y1-Y2)^2 + dz^2)
  cos_theta <- -dz / magnitude

  if (downward) cos_theta = -cos_theta

  angle_degree <- acos(pmin(pmax(cos_theta, -1), 1)) * 180/pi

  f <- function(x){ y = exp(log(100)/100*x); y[x>100]=100; y }
  net$cost <- net$cost * f(angle_degree)

  # Optional wood/leaf adjustment costs
  if (!is.null(wood_mask) & !is.null(cost_factors))
  {
    is_wood1 <- !wood_mask[net$from]
    is_wood2 <- !wood_mask[net$to]
    wood2wood <- is_wood1 & is_wood2
    leaf2leaf <- !is_wood1 & !is_wood2
    wood2leaf <- is_wood1 & !is_wood2
    net$cost[wood2wood] = net$cost[wood2wood] * cost_factors$wood2wood
    net$cost[leaf2leaf] = net$cost[leaf2leaf] * cost_factors$leaf2leaf
    net$cost[wood2leaf] = net$cost[wood2leaf] * cost_factors$wood2leaf
  }

  return(net)
}

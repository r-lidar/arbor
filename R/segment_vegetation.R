#' @export
#' @importFrom Rcpp sourceCpp
segment_vegetation = function(las, seeds, res = 0.08, k = 10, max_gap = 0.2, z_factor = 0.8)
{
  ti = tic()

  # Other parameters hard coded
  upward_cost_factor = 100
  wood2wood_cost_fact = 0.1
  leaf2leaf_cost_factor = 20
  wood2leaf_cost_factor = 1000

  # The point cloud must have hag, anisotropy and foliage computed
  attributes = names(las)
  stopifnot("anisotropy" %in% attributes)
  stopifnot("foliage" %in% attributes)
  stopifnot("foliage" %in% attributes)

  seed_coordinates <- sf::st_coordinates(seeds)
  seeds <- as.data.frame(seeds)
  seeds$X <- seed_coordinates[, "X"]
  seeds$Y <- seed_coordinates[, "Y"]
  seeds$Z <- seed_coordinates[, "Z"]

  x = mean(las$X)
  y = mean(las$Y)
  z = mean(las$Z)-100
  master_seed = data.frame(X = x, Y = y, Z = z)
  quantize(master_seed[["X"]], 0.01, las@header$`X offset`)
  quantize(master_seed[["Y"]], 0.01, las@header$`Y offset`)
  quantize(master_seed[["Z"]], 0.01, las@header$`Z offset`)
  header = rlas::header_create(master_seed)
  master_seed = suppressWarnings(lidR::LAS(master_seed, header))

  cat("Decimating the point cloud... (Step 1/6)\n")

  t0 = tic()

  dec = lidR::decimate_points(las, lidR::barycenter_per_voxel(res))

  # Calculate translation to origin
  x_translation <- mean(las$X)
  y_translation <- mean(las$Y)

  dec$X <- dec$X - x_translation
  dec$Y <- dec$Y - y_translation
  dec$Z <- dec$Z * z_factor

  header = rlas::header_create(seeds)
  seeds = lidR::LAS(seeds[, c("X", "Y", "Z", "treeID")], header)
  seeds$X <- seeds$X - x_translation
  seeds$Y <- seeds$Y - y_translation
  seeds$Z <- seeds$Z * z_factor

  num_trees <- nrow(seeds)
  num_points <- nrow(dec)

  toc(t0)

  cat("Identifying neighbors... (Step 2/6)\n") ; t0 = tic()

  # Each point is connected to its knn. The connection is bidirectional. The cost of the connection
  # is based on the euclidean distance with some variation in order to specifically follow the wood
  # and not foliage. The cost is not the same in both directions

  point_network <- compute_network(dec, k = k)
  gc()

  # Gaps have an infinite cost
  gaps = point_network$cost > max_gap
  point_network$cost[gaps] = Inf
  free(gaps)

  # It is more expensive to move in large steps
  point_network$cost <- point_network$cost^2

  # Tt is more expensive to move upward (downward actually because we are starting from the seeds)
  Z1 = dec$Z[point_network$from]
  Z2 = dec$Z[point_network$to]
  dZ = Z1-Z2
  upward = dZ > 0
  point_network$cost[upward] = point_network$cost[upward]*upward_cost_factor
  free(Z1, Z2, dZ, upward)

  is_wood1 = !dec$foliage[point_network$from]
  is_wood2 = !dec$foliage[point_network$to]
  wood2wood = is_wood1 & is_wood2
  leaf2leaf = !is_wood1 & !is_wood2
  wood2leaf = is_wood1 & !is_wood2
  point_network$cost[wood2wood] = point_network$cost[wood2wood] * wood2wood_cost_fact
  point_network$cost[leaf2leaf] = point_network$cost[leaf2leaf] * leaf2leaf_cost_factor
  point_network$cost[wood2leaf] = point_network$cost[wood2leaf] * wood2leaf_cost_factor
  free(is_wood1, is_wood2, wood2wood, wood2leaf)

  points_ids = 1:num_points

  # Each seed is connected to the knn in the point cloud. The connection is undirectional.
  # it is possible to move from the seed to the scene not the oppositite. The cost is the distance
  seed_network  <- compute_network(dec, seeds, k = k)
  seed_network$from <- seed_network$from + num_points
  seeds_ids <- (num_points+1):(num_points+num_trees)

  master_seed_network  <- compute_network(seeds, master_seed, k = num_trees)
  master_seed_network$from <- master_seed_network$from + max(seeds_ids)
  master_seed_network$to <- master_seed_network$to + min(seeds_ids)
  master_seed_network$cost <- 0.001
  master_seed_id = master_seed_network[1,1]

  toc(t0)

  cat("Constructing the graph object... (Step 3/6)\n") ; t0 = tic()

  # Extract unique tree location IDs

  combined_network <- rbind(point_network, seed_network, master_seed_network)
  free(point_network, seed_network)

  #graph_object <- cppRouting::makegraph(combined_network, directed = TRUE)
  #free(combined_network)

  toc(t0)

  from = master_seed_id
  to = points_ids

  cat("Calculating shortest paths from tree origins (can take a few minutes)... (Step 4/6)\n") ; t0 = tic()

  #distance_matrix <- cppRouting::get_distance_matrix(graph_object, from = from, to = to)
  distance_matrix <- lidRtls:::get_distance_matrix(combined_network, from, to)
  distance_matrix[is.infinite(distance_matrix)] = NA_real_

  i = colMins(distance_matrix)
  j = which(!is.na(i))

  gc()

  from = rep(master_seed_id, length(j))
  to = points_ids[j]

  t = tic()
  # For loop by chunk to reduce memory usage and have an estimated progression
  chunk_size <- 1000000
  chunks <- split(to, ceiling(seq_along(to) / chunk_size))
  pb <- utils::txtProgressBar(min = 0, max = length(chunks), style = 3, width = 50)

  treeID <- rep(NA_integer_, npoints(dec))

  for (i in seq_along(chunks))
  {
    current_to <- chunks[[i]]

    #path <- cppRouting::get_path_pair(graph_object, from = from[seq_along(current_to)], to = current_to)
    #path <- lapply(path, function(x) as.integer(x)[-1])
    #path <- unname(do.call(c, path))

    path = lidRtls:::findPaths(combined_network, from[seq_along(current_to)], current_to)
    path = path$paths
    path <- lapply(path, function(x) x[2])
    tree_id_vector <- unlist(path)
    tree_id_vector = tree_id_vector - min(seeds_ids)
    treeID[current_to] = seeds$treeID[tree_id_vector]
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  toc(t)

  free(combined_network)

  # chunk_size <- 100000
  # chunks <- split(to, ceiling(seq_along(to) / chunk_size))
  # pb <- txtProgressBar(min = 0, max = length(chunks), style = 3, width = 50)
  # treeID <- rep(NA_integer_, npoints(dec))
  #
  # for (i in seq_along(chunks))
  # {
  #   current_to <- chunks[[i]]
  #   t = tic()
  #   distance_matrix <- cppRouting::get_distance_matrix(graph_object, from = seeds_ids, to = current_to)
  #   toc(t)
  #
  #   tree_id_vector = apply(distance_matrix, 2, which.min)
  #   ids = sapply(tree_id_vector, function(x) {  if (length(x) == 0) return(NA_integer_) else return(x) })
  #   ids = unname(unlist(ids))
  #   treeID[current_to] = seeds$treeID[ids]
  #
  #   setTxtProgressBar(pb, i)
  # }
  # close(pb)

  dec <- lidR::add_lasattribute(dec, treeID, name = "treeID", desc = "tree ID")

  toc(t0)

  #plot(dec, color = "treeID", legend = TRUE)

  cat("Assigning tree IDs to the dense point cloud... (Step 6/6)\n") ; t0 = tic()

  dec$X <- dec$X + x_translation
  dec$Y <- dec$Y + y_translation
  dec$Z <- dec$Z / z_factor

  las = expand_treeid_to_neighbors(las, dec, z_factor = z_factor)

  toc(ti, space = "")

  gc()

  return(las)
}

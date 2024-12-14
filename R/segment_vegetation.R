#' @export
#' @importFrom Rcpp sourceCpp
segment_vegetation = function(las, seeds, res = 0.08, k = 10, max_gap = 0.2, z_factor = 0.8)
{
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

  toc(t0)

  cat("Constructing the graph object... (Step 3/6)\n") ; t0 = tic()

  # Extract unique tree location IDs

  combined_network <- rbind(point_network, seed_network)
  free(point_network, seed_network)

  graph_object <- cppRouting::makegraph(combined_network, directed = TRUE)
  free(combined_network)

  toc(t0)

  from = seeds_ids
  to = points_ids

  cat("Calculating shortest paths from tree origins (can take a few minutes)... (Step 4/6)\n") ; t0 = tic()

  distance_matrix <- cppRouting::get_distance_matrix(graph_object, from = from, to = to)
  gc()

  tree_id_vector = colMins(distance_matrix)
  free(distance_matrix)
  treeID <- rep(NA_integer_, npoints(dec))
  treeID[points_ids] = seeds$treeID[tree_id_vector]
  free(tree_id_vector)

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

  toc(t0, units = "mins", space = "")

  gc()

  return(las)
}

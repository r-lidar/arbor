#' @export
segment_foliage = function(las, dtm, res = 0.08, k = 10, max_gap = 0.2, min_passage = 5, th_anisotropy = 0.75)
{
  # Other parameters hard coded
  upward_cost_factor = 100
  connected_components_res = 0.05
  connected_components_min = 1000
  target_sample_res = 0.4
  wood_assignation_k = 50
  wood_assignation_dist = 0.05
  wood_extra_reasignation_k = 10
  wood_extra_reasignation_dist = 0.03

  las@data$pointID = 1:npoints(las)

  ti = tic()

  # The point cloud must have hag and anisotropy computed
  attributes = names(las)
  stopifnot("anisotropy" %in% attributes)
  stopifnot("hag" %in% attributes)

  cat("Point cloud decimation... (Step 1/9)\n") ; t0 = tic()

  dec = lidR::decimate_points(las, lidR::barycenter_per_voxel(res))

  # Calculate translation to origin for computation stability in the path finder
  x_translation <- mean(las$X)
  y_translation <- mean(las$Y)

  # The decimated points cloud of the scene
  dec$X <- dec$X - x_translation
  dec$Y <- dec$Y - y_translation
  dec$Z <- dec$Z * z_factor
  dec@header@VLR$Extra_Bytes = NULL
  dec@data = dec@data[, .(X,Y,Z, anisotropy, pointID)]
  num_points <- npoints(dec)

  # The target points spread on the volume
  target =  lidR::decimate_points(las, lidR::barycenter_per_voxel(target_sample_res))
  target$X <- target$X - x_translation
  target$Y <- target$Y - y_translation
  target$Z <- target$Z * z_factor
  target@data = target@data[, .(X,Y,Z, anisotropy, pointID)]
  num_target = npoints(target)

  # A layer of ground points used as a connector to the unique seed
  gnd = seed_from_dtm(dtm)
  gnd$X <- gnd$X - x_translation
  gnd$Y <- gnd$Y - y_translation
  gnd$Z <- gnd$Z * z_factor
  quantize(gnd[["X"]], 0.01, las@header$`X offset`)
  quantize(gnd[["Y"]], 0.01, las@header$`Y offset`)
  quantize(gnd[["Z"]], 0.01, las@header$`Z offset`)
  gnd$anisotropy = 1
  gnd$pointID = 0
  header <- rlas::header_create(gnd)
  gnd = suppressWarnings(lidR::LAS(gnd, header))
  num_gnd <- npoints(gnd)

  # The seed
  x = mean(gnd$X)
  y = mean(gnd$Y)
  z = mean(gnd$Z)-100
  seed = data.frame(X= x, Y = y, Z = z, anisotropy = 1, pointID = 0)
  quantize(seed[["X"]], 0.01, las@header$`X offset`)
  quantize(seed[["Y"]], 0.01, las@header$`Y offset`)
  quantize(seed[["Z"]], 0.01, las@header$`Z offset`)
  header = rlas::header_create(seed)
  seed = suppressWarnings(lidR::LAS(seed, header))

  # Bind the cloud with the target point for the path finder
  #point_coordinates = rbind(dec, target)

  toc(t0)

  cat("Building point cloud connectivity... (Step 2/9)\n") ; t0 = tic()

  # Each point is connected to its knn. The connection is bidirectional. The cost of the connection
  # is based on the euclidean distance with some variation in order to specifically follow the wood
  # and not foliage. The cost is not the same in both directions

  point_network <- compute_network(dec, k = k)
  gc()

  # Gaps have an infinite cost
  gaps = point_network$cost > max_gap
  point_network$cost[gaps] = Inf

  # It is more expensive to move in large steps
  point_network$cost <- point_network$cost^2

  # It is more expensive to move upward  (downward actually because we are starting from the seeds)
  Z1 = dec$Z[point_network$from]
  Z2 = dec$Z[point_network$to]
  upward = Z1-Z2 > 0
  point_network$cost[upward] = point_network$cost[upward]*upward_cost_factor
  free(upward, Z1, Z2)

  # The cost is weighted by the anisotropy
  #A1 = point_coordinates$anisotropy[point_network$from]
  #A2 = point_coordinates$anisotropy[point_network$to]
  #W = 1-((A1+A2)/2)
  #W <- 0.1 + W * (0.9 - 0.1)
  #point_network$cost = point_network$cost * W
  #free(W, A1, A2)

  points_ids = 1:num_points

  toc(t0)

  cat("Building target connectivity... (Step 3/9)\n")

  # Point cloud is connected to targets to reach. The connection is undirectional. It is possible to
  # move from the point cloud to the targets not the opposite. The cost of the connection is
  # virtually null because the target are subsampled from the point cloud

  target_network <- compute_network(dec, target, k = 1)
  from = target_network$from
  target_network$from = target_network$to # switch direction
  target_network$to = from
  target_network$to <- target_network$to + num_points
  target_ids = (num_points+1):(num_points+1+num_target)

  cat("Building ground connectivity... (Step 3/9)\n") ; t0 = tic()

  # The ground points are connected to their knn points in the scene. The connection is undirectional.
  # it is possible to move from the ground to the scene not the oppositite. The cost of the connection
  # is the euclidean distance.

  ground_network = compute_network(dec, gnd, k = k)
  ground_network$from <- ground_network$from + num_points + num_target
  ground_ids = (num_points+num_target+1):(num_points+num_target+1+num_gnd)

  toc(t0)

  cat("Building connectivity of a single seed to the ground... (Step 4/9)\n") ; t0 = tic()

  # The unique seed is connected to all the ground points in the scene. The connection cost is constant
  # and virtually 0. It allows to start from a single point to reach all the targets

  seed_network  <- compute_network(gnd, seed, k = num_gnd)
  seed_network$from <- seed_network$from + max(ground_ids)
  seed_network$to <- seed_network$to + min(ground_ids)
  seed_network$cost <- 0.001
  seed_id = seed_network[1,1]

  toc(t0)

  cat("Constructing the graph object... (Step 5/9)\n") ; t0 = tic()

  combined_network <- rbind(point_network, target_network, ground_network, seed_network)
  free(point_network, target_network, ground_network, seed_network)

  #graph_object <- cppRouting::makegraph(combined_network, directed = TRUE)

  #free(combined_network)

  toc(t0)

  cat("Calculating shortest paths to ground... (Step 6/9)\n") ; t0 = tic()

  # Shortest distance in order to find reachable 'target' points
  # I don't remember why we need that. In theory cppRouting::get_path_pair() should be enough but maybe non reachable targets pose an issue
  #distance_matrix <- cppRouting::get_distance_matrix(graph_object, from = seed_id, to = target_ids, allcores = FALSE)
  distance_matrix <- lidRtls:::get_distance_matrix(combined_network, seed_id, target_ids)
  distance_matrix[is.infinite(distance_matrix)] = NA_real_

  i = colMins(distance_matrix)
  j = which(!is.na(i))

  free(distance_matrix, i)
  toc(t0)

  cat("Calculating paths (can take a few minutes)... (Step 7/9)\n") ; t0 = tic()

  dec@data$count = 0

  from = rep(seed_id, length(j))
  to = target_ids[j]

  # For loop by chunk to reduce memory usage and have an estimated progression
  chunk_size <- 50000
  chunks <- split(to, ceiling(seq_along(to) / chunk_size))
  pb <- utils::txtProgressBar(min = 0, max = length(chunks), style = 3, width = 50)

  for (i in seq_along(chunks))
  {
    current_to <- chunks[[i]]

    #path <- cppRouting::get_path_pair(graph_object, from = from[seq_along(current_to)], to = current_to)
    #path <- lapply(path, function(x) as.integer(x)[-1])
    #path <- unname(do.call(c, path))

    path = lidRtls:::findPaths(combined_network, from[seq_along(current_to)], current_to)
    path = path$paths
    path <- lapply(path, function(x) x[-1])
    path <- do.call(c, path)

    count <- table(path)
    id <- as.numeric(names(count))

    rm <- id > num_points
    id <- id[!rm]
    count <- count[!rm]

    dec@data$count[id] <- dec@data$count[id] + count

    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  free(combined_network)


  dec@data$count[dec@data$count > 0] = log(dec@data$count[dec@data$count > 0])
  #x = plot(dec, color = "count", legend = T)
  #plot(dec, color = "anisotropy", legend = T, breaks = "quantile")

  free(path, count, id, rm, from, to)

  toc(t0)

  cat("Assigning wood to small structure... (Step 8/9)\n") ; t0 = tic()

  skeleton = filter_poi(dec, count > log(min_passage))

  # For visualization and debugging mainly
  las@data$skeleton = 0
  las@data$skeleton[skeleton$pointID] = 1
  las = add_lasattribute_manual(las, name = "skeleton", desc = "skeleton points", type = "char")
  #plot(las, color = "skeleton")

  #x = plot(dec, color = "count", legend = T)
  #plot(skeleton, add = x, pal = "red", size = 2)


  # The decimated points
  skeleton$X <- skeleton$X + x_translation
  skeleton$Y <- skeleton$Y + y_translation
  skeleton$Z <- skeleton$Z / z_factor

  skeleton_neighbors  <- knnx(las, skeleton, k = wood_assignation_k)
  rm = skeleton_neighbors$nn.dist > wood_assignation_dist
  id = skeleton_neighbors$nn.index[!rm]

  las@data$wood = FALSE
  las@data$wood[id] = TRUE

  free(skeleton_neighbors, rm, id)

  #plot(las, color = "wood", pal =c("gray", "chocolate4"))

  toc(t0)

  #x = plot(las)
  #plot(skeleton, pal = "red", size = 4, add = x)

  cat("Filter anisotropy... (Step 7/6)\n") ; t0 = tic()

  nofoliage = filter_poi(las, anisotropy > th_anisotropy | wood == TRUE)
  #foliage = filter_poi(las, anisotropy <= th_anisotropy & wood == FALSE)

  toc(t0)

  cat("Connected component cleaning... (Step 8/9)\n") ; t0 = tic()

  nofoliage = connected_components(nofoliage, connected_components_res, connected_components_min)
  #foliage2 = nofoliage[nofoliage$clusterID == 0 & nofoliage$wood == FALSE]
  nofoliage = nofoliage[nofoliage$clusterID != 0 | nofoliage$wood == TRUE]

  toc(t0)

  #x = plot(nofoliage, pal = "chocolate4")
  #plot(foliage, pal = "forestgreen", add = x)
  #plot(foliage2, pal = "forestgreen", add = x)

  cat("Extra wood reasignation... (Step 9/9)\n") ; t0 = tic()

  wood_neighbors  <- knnx(las, nofoliage, k = wood_extra_reasignation_k)
  rm = wood_neighbors$nn.dist > wood_extra_reasignation_dist
  id = wood_neighbors$nn.index[!rm]

  las@data$foliage = TRUE
  las@data$foliage[id] = FALSE
  las@data$foliage = as.integer(las$foliage)
  las = add_lasattribute_manual(las, name = "foliage", desc = "foliage 1 wood 0", type = "char")

  free(nofoliage, rm, id, wood_neighbors)

  toc(t0)

  toc(ti, space = "")

  gc()

  return(las)
}

#' Find seeds to perform instance segmentation
#'
#' In order to perform instance segmentation with  \link{segment_instance} we need
#' some seeds with reference treeIDs. Thus function finds the seeds
#'
#' @param las A LAS object from lidR.
#' @param params list See \link{parameters}.
#' @export
find_seeds <- function(las, params)
{
  ti = tic()

  foliage <- clusterID <- max_diameter <- passage <- hag <- X <- Y <- Z <- . <- NULL

  # The point cloud is supposed to have passage hag and foliage
  attributes <- names(las)
  stopifnot("hag" %in% attributes)
  stopifnot("passage" %in% attributes)
  stopifnot("foliage" %in% attributes)

  # The heights at which we slice
  # Extract some slices of wood (thickness 3cm)
  thick  <- params$seed$slice_thickness
  heights <- c(min(las$hag)+thick, params$seed$slice_at)

  # Extract the passages of big trees and small features
  th                   <- max(heights) + 0.2
  min_passage          <- params$seed$min_passage
  long_passages        <- lidR::filter_poi(las, passage > min_passage, hag < th)
  short_passages       <- lidR::filter_poi(las, passage > 1, hag < min(hag) + 0.5, passage < 10)
  long_passages@data   <- long_passages@data[, .(X,Y,Z, passage, hag)]
  short_passages@data  <- short_passages@data[, .(X,Y,Z, passage, hag)]
  if (FALSE)
  {
    plot(long_passages)
    plot(short_passages)
  }

  cat("Slicing the point cloud...\n") ; t0 = tic()

  slices <- slice_poi(las, heights, thick)
  wood   <- lidR::filter_poi(slices, foliage == 0)

  # wood can be noisy we need to clean noise.
  #sor_k     <- params$seed$sor_k
  #sor_m     <- params$seed$sor_m
  #wood      <- lidR::classify_noise(wood, lidR::sor(k = sor_k, m = sor_m)) # very very aggressive sor
  #wood      <- lidR::remove_noise(wood)
  wood@data <- wood@data[, .(X,Y,Z, passage, hag)]

  if (FALSE)
  {
    x = plot(wood)
    plot(long_passages, add = x, pal = "green")
  }

  toc(t0, "  ")
  cat("Connected component...\n") ; t0 = tic()

  # Connect the wood point into clusters
  cl_wood <- connected_components(wood, 0.05, 10, connectivity = 26)
  cl_wood <- lidR::filter_poi(cl_wood, clusterID != 0)

  toc(t0, "  ")
  cat("Circle detection per component...\n") ; t0 = tic()

  # For each cluster search for circles. If we have a nice circle we have a tree
  is.valid.circle <- function(radius, angle_range, pinliner, pinside)
  {
    if (radius < 0.02)  return(FALSE)
    if (radius  < 0.05)  return(angle_range > 180 & pinliner > 30)
    if (pinside > 20)   return(FALSE)
    if (radius < 0.10)  return(angle_range > 90 & pinliner > 50)
    return(angle_range > 130 & pinliner > 30)
  }
  fit_circle_to_seed <- function(cl)
  {
    id = cl$clusterID[1]
    if (nrow(cl) < 20) return(NULL)
    cl <- as.matrix(cl[,1:3])
    circle <- ransac_circle(cl, num_iterations = 400, inlier_threshold = 0.02)

    valid  <- is.valid.circle(circle$radius, circle$covered_arc_degree, circle$percentage_inlier*100, circle$percentage_inside*100)

    if (FALSE)
    {
      if (valid) col = "darkgreen" else col = "red"
      plot(cl$X, cl$Y, asp = 1, main = id)
      symbols(circle$center_x, circle$center_y, circles = circle$radius, inches = FALSE, add = TRUE, fg = col)
      symbols(circle$center_x, circle$center_y, circles = circle$radius+0.01, inches = FALSE, add = TRUE, fg = col)
      symbols(circle$center_x, circle$center_y, circles = circle$radius-0.01, inches = FALSE, add = TRUE, fg = col)
    }

    if (valid) return(data.frame(X = circle$center_x, Y = circle$center_y, Z = circle$z, R = circle$radius, id = id))
    else return(NULL)
  }

  clusters <- split(cl_wood@data, by = "clusterID")
  clusters <- Filter(function(x) { nrow(x) > 20 }, clusters)

  n <- length(clusters)
  i <- 0
  circles <- lapply(clusters, function(cl)
  {
    i <<- i + 1
    if (i %% 10 == 0)  cat(sprintf("\r  Processed %d / %d", i, n))
    fit_circle_to_seed(cl)
  })

  cat("\n")
  circles <- Filter(Negate(is.null), circles)

  circles_detected = length(circles) > 0

  if (!circles_detected)
    warning("No circle dectected")
  else
    circles  <- do.call(rbind, circles)

  toc(t0, "  ")

  if (FALSE)
  {
    x <- plot(cl_wood, color = 'clusterID', pal = pastel.colors(500))
    plot(long_passages, pal = "green", add = x)
    for (i in 1:nrow(circles))
      add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], circles$Z[i])
  }

  cat("Safe zone exclusion...\n") ; t0 = tic()

  if (circles_detected)
  {
    # For each circle we exclude wood in a safe zone beyond the circles.
    # This allow to clean false positive around important trees and prevent
    # dummy connection caused by noise in the next connected component step
    px <- wood$X
    py <- wood$Y
    pz <- wood$Z
    rm <- rep(FALSE, lidR::npoints(wood))
    safe_zone <- 0.2
    for (i in 1:nrow(circles))
    {
      cx <- circles$X[i]
      cy <- circles$Y[i]
      cz <- circles$Z[i]
      r  <- circles$R[i]
      d  <- sqrt((px-cx)^2 + (py-cy)^2 +(pz-cz)^2)
      rm[d > (r + 0.02) & d  < (r + safe_zone)] <- TRUE
    }
    wood <- wood[!rm]

    toc(t0, "  ")

    if (FALSE)
    {
      plot(wood)
      for (i in 1:nrow(circles))
        add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], circles$Z[i])
    }

    if (FALSE)
    {
      x = plot(long_passages, pal = "green")
      plot(short_passages, add = x, pal = "red")
      plot(wood, add = x, pal = foliage.colors[1])
      for (i in 1:nrow(circles))
        add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], circles$Z[i])
    }

    cat("Overlapping disc detection...\n") ; t0 = tic()

    # Overlapping discs
    # ------------------
    # Pairwise distances between centers
    circles$id = NULL
    df = circles
    dist_mat <- as.matrix(stats::dist(df[, c("X", "Y")], diag = TRUE, upper = TRUE))

    # Define an overlap rule (e.g. centers closer than sum of radii)
    overlap <- dist_mat < outer(df$R, df$R, "+") & dist_mat > 0

    # Build adjacency graph from overlap matrix
    n <- nrow(df)
    visited <- logical(n)
    group <- integer(n)
    gid <- 0

    for (i in seq_len(n)) {
      if (!visited[i]) {
        gid <- gid + 1
        # breadth-first search (BFS) for connected components
        queue <- i
        while (length(queue) > 0) {
          node <- queue[1]
          queue <- queue[-1]
          if (!visited[node]) {
            visited[node] <- TRUE
            group[node] <- gid
            neighbors <- which(overlap[node, ])
            queue <- c(queue, neighbors[!visited[neighbors]])
          }
        }
      }
    }

    df$group_id <- group
    df
    circles = df

    toc(t0, "  ")
    cat("Discs to seeds conversion...\n") ; t0 = tic()

    # Convert circle to points
    res= params$decimation$barycentric_predecimation_resolution
    circle_points_list <- lapply(1:nrow(circles), function(i) {
      generate_circle_points(circles$X[i], circles$Y[i], circles$Z[i], circles$R[i], step = res)
    })
    circle_points = do.call(rbind, circle_points_list)

    connectors <- generate_all_groups(circles, step_z = res)
    connectors = rbind(connectors$disks)

    circle_points = rbind(circle_points, connectors)
    circle_points$passage = 1000
    circle_points$hag = 0
    circle_points = suppressWarnings(lidR::LAS(circle_points, lidR::header(las)))

    # Bind the wood and the long passages and compute connected component and merge passage from the same trees
    lidR::st_crs(long_passages) = lidR::st_crs(wood)
    lidR::st_crs(circle_points) = lidR::st_crs(wood)
    temp   <- suppressWarnings(rbind(wood, long_passages, circle_points))
    res    <- round(params$decimation$barycentric_predecimation_resolution*0.8, 2)
    temp$Z <- temp$Z * 0.5
    temp   <- connected_components(temp, res, 1, name = "treeID", connectivity = 26)
    temp$Z <- temp$Z / 0.5

    if (FALSE)
    {
      x <- plot(temp, color = "treeID")
      #plot(passage, add = x)
      for (i in 1:nrow(circles))
        add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], circles$Z[i])
    }
  }

  toc(t0, "  ")
  cat("Pathfinder...\n") ; t0 = tic()

  # We have seed for the big trees (long passage)
  long_passages_seeds = lidR::filter_poi(temp, passage > 0)

  if (FALSE)
  {
    x <- plot(long_passages_seeds, color = "treeID")
    plot(short_passages, add = x, pal = "red")
  }

  # Some short passage could be from big tree. Path finder to attach them.
  short_passages@data$foliage = 0
  p = default_arbor_parameters
  p$path_finder$max_gap = 0.1
  p$path_finder$k_neighborhood_connectivity = 10
  p$path_finder$k_seed_connectivity = 2
  p$path_finder$distance_power = 1
  p$path_finder$angle_penalty = function(x) {return(rep(1, length(x))) }

  sink(tempfile())
  on.exit(suppressWarnings(sink()), add = TRUE)
  short_passages = segment_vegetation(short_passages, long_passages_seeds, p)
  sink()

  if (FALSE)
  {
    u = short_passages
    u@data$foliage = NULL
    u@data$pointID = NULL
    v = rbind(long_passages_seeds, u)
    plot(v, color = "treeID")
  }

  toc(t0, "  ")
  cat("Connected component...\n") ; t0 = tic()

  # Some short passage don't have IDs
  short_passages_noid = short_passages[is.na(short_passages$treeID)]
  short_passages_noid = lidR::connected_components(short_passages_noid, 0.1, 1, "treeID")
  short_passages_noid$treeID = short_passages_noid$treeID + max(long_passages_seeds$treeID, na.rm = TRUE) + 1


  # Bind the wood and the passage in order to compute
  # connected component and merge passage from the same trees for big tree only
  short_passages_withid = short_passages[!is.na(short_passages$treeID)]
  short_passages_withid$foliage = NULL
  short_passages_withid$pointID = NULL
  short_passages_noid$foliage = NULL
  short_passages_noid$pointID = NULL

  toc(t0, "  ")

  if (FALSE)
  {
    x = plot(rbind(long_passages_seeds, short_passages_withid), color = "treeID")
    plot(short_passages_noid, add = x, pal = "red")
  }

  seeds = suppressWarnings(rbind(long_passages_seeds, short_passages_withid, short_passages_noid))

  # Retain only the seed below BH
  seeds <- lidR::filter_poi(seeds, hag < 1)

  toc(ti, "")

  seeds
}

slice_poi = function(las, heights, thinkness = 0.02)
{
  # Build dynamic filter for slices
  slice_filter <- Reduce(`|`, lapply(heights, function(s)
  {
    (las$hag > (s-thinkness/2) & las$hag < (s + thinkness/2))
  }))

  lidR::filter_poi(las, slice_filter)
}

# Convert to points
generate_circle_points <- function(x, y, z, r, step = 0.1)
{
  # Circumference
  circumference <- 2 * pi * r
  # Number of points for approximately every `step` meters
  n_points <- ceiling(circumference / step)
  # Angles for points
  theta <- seq(0, 2*pi, length.out = n_points + 1)[-1]  # remove last point to avoid duplicate
  # Generate points
  data.frame(
    X = x + r * cos(theta),
    Y = y + r * sin(theta),
    Z = rep(z, n_points)
  )
}

# --- Generate random points on a disk (XY plane) ---
generate_disk_points <- function(x, y, z, r, n = 1000)
{
  theta  <- stats::runif(n, 0, 2*pi)
  radius <- sqrt(stats::runif(n)) * r
  data.frame(
    X = x + radius * cos(theta),
    Y = y + radius * sin(theta),
    Z = rep(z, n)
  )
}

generate_disk_radii <- function(x, y, z, r, n = 8)
{
  angles <- seq(0, 2*pi, length.out = n + 1)[- (n + 1)]  # remove last to avoid duplicating 0

  data.frame(
    X = x + r * cos(angles),
    Y = y + r * sin(angles),
    Z = rep(z, n)
  )
}


# --- Generate points for one group ---
generate_group_points <- function(df_group, step_z = 0.05)
{
  df_group <- df_group[order(df_group$Z), ]
  all_points <- list()
  centerline <- list()

  for (i in seq_len(nrow(df_group) - 1)) {
    c1 <- df_group[i, ]
    c2 <- df_group[i + 1, ]

    # vertical interpolation sequence (every step_z)
    z_seq <- seq(c1$Z, c2$Z, by = step_z)
    if (utils::tail(z_seq, 1) != c2$Z)
      z_seq <- c(z_seq, c2$Z)  # include top

    t_seq <- (z_seq - c1$Z) / (c2$Z - c1$Z)

    # interpolate center and radius
    x_seq <- c1$X + t_seq * (c2$X - c1$X)
    y_seq <- c1$Y + t_seq * (c2$Y - c1$Y)
    r_seq <- c1$R + t_seq * (c2$R - c1$R)

    # add centerline points
    centerline[[length(centerline) + 1]] <- data.frame(X = x_seq, Y = y_seq, Z = z_seq)

    # add 4-radii disks
    for (j in seq_along(z_seq)) {
      all_points[[length(all_points) + 1]] <- generate_disk_radii(
        x_seq[j], y_seq[j], z_seq[j], r_seq[j]
      )
    }
  }

  list(
    disks = do.call(rbind, all_points),
    centerline = do.call(rbind, centerline)
  )
}

# --- Wrapper for all groups ---
generate_all_groups <- function(df, step_z = 0.05) {
  groups <- unique(df$group_id)
  res_disks <- list()
  res_centers <- list()

  for (g in seq_along(groups)) {
    df_g <- df[df$group_id == groups[g], ]
    r <- generate_group_points(df_g, step_z)
    res_disks[[g]] <- r$disks
    res_centers[[g]] <- r$centerline
  }

  list(
    disks = do.call(rbind, res_disks),
    centerline = do.call(rbind, res_centers)
  )
}

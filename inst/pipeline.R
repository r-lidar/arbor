library(lidR)
library(lidRtls)

z_factor = 0.8
k = 10
k_ani = 75
max_gap = 0.2
ncpu = 12
slice_seeds_at = c(0.5, 1)
path_finder_resolution = 0.1
voxel_sampling_resolution = 0.01
cut_above_ground = 0.5

display = FALSE

select = "tircn0"
filter = ""
filter = "-keep_random_fraction 0.2"
file = "~/Documents/Entreprise/clients/fsinvestor/SanDiego/FSTESTSCAN3UCSD_01_laz1_4_extract30m.laz"
file = "~/Documents/Entreprise/clients/fsinvestor/ST-X-ZamPlot/ZamPlot_part1.laz"
file = "~/Documents/Entreprise/clients/fsinvestor/ST-X-ZamPlot/ZamPlot_part2.laz"
file = "~/Documents/Entreprise/clients/fsinvestor/ST-X-ZamPlot/ZamPlot_part3.laz"
file = "~/Documents/Entreprise/clients/fsinvestor/Rwanda/Kwandahillside/Kwandahillside.laz"
file = "~/Documents/Entreprise/clients/fsinvestor/Rwanda/Referencesite1_part1.laz"
file = "~/Documents/Entreprise/clients/fsinvestor/Rwanda/Referencesite1_part2.laz"



# ====== READ POINT CLOUD =======

las = readTLS(file, select = select, filter = filter)
las

#las = clip_circle(las, 30,25,5)
#las = clip_circle(las, 70,-50,5)

# ====== REDUCE TO 1 POINTS/CM3 =======

algo = lowest_attribute_per_voxel(voxel_sampling_resolution, "Range")
las = decimate_points(las, algo)
if (display) plot(las)

#las = lidR::knn_distance(las, k = k)
#plot(las, color = "distance", breaks = "quantile", legend = T)


# ====== GROUND CLASSIFICATION ====

las = classify_ground(las, lidR::csf(rigidness = 1, class_threshold = 0.05, cloth_resolution = 0.1), last_returns = FALSE)
ground = filter_poi(las, Classification == LASGROUND)
ground = decimate_points(ground, lowest(0.25))
ground = classify_noise(ground, sor(k = 10, m = 2))
ground = remove_noise(ground)
ground$Classification = lidR::LASGROUND

# ====== DTM & HAG =======

dtm = rasterize_terrain(ground, 0.5, tin())
las = height_above_ground(las, algorithm = tin(), dtm = dtm)

# ====== KEEP ABOVE DTM ======

bottom = filter_poi(las, hag <= cut_above_ground)
las = filter_poi(las, hag > cut_above_ground)

if (display) plot(las)

# ===== COMPUTE ANISOTROPY =======

olas = sf::st_coordinates(las)

las = lidRtls:::smooth3d(las, 0.04)
las = lidR::knn_distance(las, k = 20)
f <- ecdf(las$distance)
las@data$anisotropy <- 1-f(las$distance)

las = compute_anisotropy(las, k = k_ani)

if (display) plot(las, color = "anisotropy", legend = T, breaks = "quantile")

# ====== SEGMENT FOLIAGE/WOOD ======

las = segment_foliage(las, dtm, res = 0.08, max_gap = max_gap, k = k, min_passage = 5, th_anisotropy = 0.75)
#las = segment_foliage(las, dtm, res = 0.02, max_gap = 0.05, k = k, min_passage = 5, th_anisotropy = 0.6, space_res = 0.1)

if (display)
{
plot(las, color = "foliage", pal = c("#523b0a", "darkgreen"), size = 2) |> add_dtm3d(dtm)
plot(filter_poi(las, foliage == FALSE), pal = "#523b0a", size = 2) |> add_dtm3d(dtm)
x = plot(las, color = "foliage", pal = c("#523b0a", "darkgreen"), size = 2) |> add_dtm3d(dtm)
sk = filter_poi(las, skeleton == T)
plot(sk, add = x, pal = "red", size = 3)
}


#writeLAS(las, "~/Téléchargements/segmented.las")
#las = readLAS("~/Téléchargements/segmented.las", select = "xyz01234")

# ====== FIND TREE SEEDS =======

seeds = find_seeds(las, slice_seeds_at)

if (display)
{
col = pastel.colors(length(unique(seeds$treeID)))
col = col[as.integer(as.factor(seeds$treeID))]
plot(filter_poi(las, hag < 3), color = "foliage", pal = c("chocolate4", "darkgreen")) |> add_treetops3d(seeds, radius = 0.1, col = col)  |> add_dtm3d(dtm)
plot(filter_poi(las, hag < slice_seeds_at[2]), color = "foliage", pal = c("chocolate4", "darkgreen")) |> add_treetops3d(seeds, radius = 0.05, color = col) |> add_dtm3d(dtm)
}

# ====== SEGMENT TREES =======

las = segment_vegetation(las, seeds, max_gap = max_gap, k = k)

if (display) x = plot(las, color = "treeID") |> add_dtm3d(dtm) |> add_treetops3d(seeds, radius = 0.08)

trees = clean_small_cluster(las, max_heigh = 5)

if (display)
{
plot(trees, color = "skeleton", legend = TRUE) |> add_dtm3d(dtm)
plot(trees, color = "treeID", legend = TRUE, size = 2) |> add_dtm3d(dtm)
plot(trees, color = "foliage", pal = c("chocolate4", "darkgreen")) |> add_dtm3d(dtm)
plot(filter_poi(trees, foliage == FALSE), color = "treeID", legend = TRUE) |> add_dtm3d(dtm) |> add_treetops3d(seeds, radius = 0.1)
}

# ==== CLEAN BOTTOM ADN TREE EXTENSION =====

las = fix_splited_trees(las)

generate_cylinder_points <- function(circle, height = 0.5, n_points = 1000)
{
  # Extract circle parameters
  center_x <- circle$center_x
  center_y <- circle$center_y
  radius <- circle$radius
  z_top <- circle$z
  z_bottom <- z_top - height

  # Create a grid of points
  cylinder_points <- data.frame(
    X = 0,
    Y = 0,
    theta = runif(n_points, 0, 2*pi),
    Z = runif(n_points, z_bottom, z_top)
  )

  # Convert to Cartesian coordinates
  cylinder_points$X <- center_x + radius * cos(cylinder_points$theta)
  cylinder_points$Y <- center_y + radius * sin(cylinder_points$theta)

  # Drop theta column (optional)
  cylinder_points$theta <- NULL

  # Return the result as a data frame
  return(cylinder_points)
}

align_to_z <- function(main_axis)
{
    # Normalize the main axis vector
    main_axis <- main_axis / sqrt(sum(main_axis^2))

    # Z-axis unit vector
    z_axis <- c(0, 0, 1)

    # Cross product to find the rotation axis
    rotation_axis <- c(
      main_axis[2] * z_axis[3] - main_axis[3] * z_axis[2],
      main_axis[3] * z_axis[1] - main_axis[1] * z_axis[3],
      main_axis[1] * z_axis[2] - main_axis[2] * z_axis[1]
    )

    # Normalize the rotation axis
    axis_length <- sqrt(sum(rotation_axis^2))
    if (axis_length > 1e-6) { # Avoid division by zero
      rotation_axis <- rotation_axis / axis_length
    } else {
      # If the main axis is already aligned with Z, return identity matrix
      return(diag(3))
    }

    # Compute the angle between main_axis and Z-axis
    angle <- acos(sum(main_axis * z_axis))

    # Construct the rotation matrix using Rodrigues' rotation formula
    K <- matrix(c(
      0, -rotation_axis[3], rotation_axis[2],
      rotation_axis[3], 0, -rotation_axis[1],
      -rotation_axis[2], rotation_axis[1], 0
    ), nrow = 3, byrow = TRUE)

    R <- diag(3) + sin(angle) * K + (1 - cos(angle)) * (K %*% K)

    return(R)
}



# bad id = 137
# bad id = 598
# bad id = 334
# bad id = 625
#trees@data$pointID = 1:npoints(trees)
extensions = list()
for (id in unique(trees$treeID))
{
  print(id)
  tt = filter_poi(trees, treeID == id, foliage == FALSE, hag < 3)
  if (is.empty(tt)) next
  tt$Z = tt$Z * 0.1
  tt = connected_components(tt, 0.05, 200)
  tt$Z = tt$Z * 10

  ids = 1

  if (length(unique(tt$clusterID)) > 1)
  {
    ids = table(tt$clusterID)
    ids = as.numeric(names(ids[which.max(ids)]))
    pid = tt$pointID[tt$clusterID != ids]
    las$foliage[pid] = TRUE
    tt = filter_poi(tt, clusterID == ids)
    #plot(tt, color = "clusterID")
    #cat("  ", id, "\n")
    #plot(tt, color = "foliage", pal = c("chocolate4", "darkgreen"))
  }

  tt = filter_poi(tt, hag < 2)

  #xyz = sf::st_coordinates(tt)
  #pca <- prcomp(xyz, center = TRUE, scale. = FALSE)
  #main_axis <- pca$rotation[, 1]  # First principal component

  # Compute the rotation matrix
  #rotation_matrix <- align_to_z(main_axis)

  # Apply the rotation to the point cloud
  #rotated_xyz = xyz %*% t(rotation_matrix)
  #rotated_xyz <- scale(rotated_xyz, center = TRUE, scale = FALSE)
  #ttt <- as.data.frame(rotated_xyz)
  #names(ttt) = c("X", "Y", "Z")
  #ttt$clusterID = tt$clusterID
  #ttt = LAS(ttt)


  # Plot the point cloud
  #rgl::plot3d(centered_xyz, col = "blue", size = 2)
  #rgl::points3d(rotated_xyz, col = "red", size = 2)
  ##rgl::arrow3d(p0 = c(0,0,0), p1 =  main_axis, type = "lines",  col = "red", length = 2)

  ranges = expand.grid(bottom = seq(0.1,1,0.05), top = seq(0.1,1,0.05))
  ranges = ranges[ranges$top - ranges$bottom >= 0.1,]
  ranges = ranges[ranges$top - ranges$bottom <= 0.6,]

  circles = apply(ranges, 1, function(x)
  {
    bottom = filter_poi(tt, hag  >= x[1], hag <= x[2])
    if (is.empty(bottom)) return(NULL)

    bottom$Z = bottom$Z * 0.01
    bottom = connected_components(bottom, 0.01, 5)
    ids = table(bottom$clusterID)
    ids = as.numeric(names(ids[which.max(ids)]))
    bottom = bottom[bottom$clusterID == ids]

    if (npoints(bottom) < 10) return(NULL)

    bottom$Z = bottom$Z * 100
    circle = fit_circle(bottom)
    inliner = length(circle$inliers)/npoints(bottom) * 100

    if (display)
    {
    plot(sf::st_coordinates(bottom), asp = 1, main = id)
    symbols(circle$center_x, circle$center_y, circles = circle$radius,  add = TRUE, fg = "red", inches = FALSE)
    symbols(circle$center_x, circle$center_y, circles = circle$radius+0.01,  add = TRUE, fg = "red", lty=3, inches = FALSE)
    symbols(circle$center_x, circle$center_y, circles = circle$radius-0.01,  add = TRUE, fg = "red", lty= 3, inches = FALSE)
    mtext(paste0("Radius = ", round(circle$radius,2), " inliner = ", round(inliner), "% sector ", circle$angle_range, " deg"))
    }
    circle$pinlier = inliner
    circle$inliers = NULL
    as.data.frame(circle)
  })

  if (is.null(circles)) next

  circles = do.call(rbind, circles)

  valid = circles$radius < 0.25 & circles$angle_range > 90
  n = sum(valid)
  if (n == 0) valid = circles$radius < 0.25
  n = sum(valid)
  if (n == 0) next

  circles = circles[valid,]
  circles = circles[rev(order(circles$pinlier)),]
  n = min(c(5, n))

  if (n == 0) next

  circles = circles[1:n,]
  x = median(circles$center_x)
  y = median(circles$center_y)
  z = min(circles$z)+0.05
  r = median((circles$radius))
  circle = list(center_x = x, center_y = y, z = z, radius = r)

  loc = matrix(c(x,y), ncol = 2)
  zgnd = terra::extract(dtm, loc, method = "bilinear")
  hcyl = as.numeric(z-zgnd+0.15)

  if (circle$radius > 0.1) cat("Big tree", id, "\n")

  if (circle$radius < 0.15)
  {
    extension = generate_cylinder_points(circle, height = hcyl)
    #extension = as.matrix(extension)
    #extension <- extension %*% rotation_matrix
    #extension = as.data.frame(extension)
    #names(extension) = c("X", 'Y', "Z")
    quantize(extension$X, tt@header[["X scale factor"]], tt@header[["X offset"]])
    quantize(extension$Y, tt@header[["Y scale factor"]], tt@header[["Y offset"]])
    quantize(extension$Z, tt@header[["Z scale factor"]], tt@header[["Z offset"]])
    extension$treeID = id
    extensions[[as.character(id)]] = extension
  }
  #x = plot(tt)
  #plot(LAS(extension), add = x)
}

combine_with_fill <- function(df1, df2, fill_value = 0L)
{
  all_cols <- union(names(df1), names(df2))

  # Add missing columns with the fill value
  for (col in setdiff(all_cols, names(df1))) {
    df1[[col]] <- fill_value
  }
  for (col in setdiff(all_cols, names(df2))) {
    df2[[col]] <- fill_value
  }

  # Ensure column order matches
  df1 <- df1[, ..all_cols]
  df2 <- df2[, ..all_cols]

  # Combine rows
  rbind(df1, df2)
}


trees = clean_small_cluster(las, max_heigh = 5)

extensions = do.call(rbind, extensions)
extensions$anisotropy = 1
extensions$wood = TRUE
data.table::setDT(extensions)
x = plot(filter_poi(trees, foliage == FALSE)) |> add_dtm3d(dtm)
plot(LAS(extensions), add = x)


trees@data = combine_with_fill(trees@data, extensions)
trees = las_update(trees)

plot(trees, color = "treeID") |> add_dtm3d(dtm)
plot(trees, color = "foliage", pal = c("chocolate4", "darkgreen")) |> add_dtm3d(dtm)



# ==== CLIP BUFFER ======

valid_trees = clip_buffer(las, seeds)
valid_trees = clean_small_cluster(valid_trees, max_heigh = 6)
plot(filter_poi(valid_trees, foliage == FALSE), color = "treeID", legend = TRUE) |> add_dtm3d(dtm) |> add_treetops3d(seeds, radius = 0.1)


o =  tools::file_path_sans_ext(file)
r = paste0(o, "_dtm.tif")
s = paste0(o, "_seeds.shp")
o = paste0(o, "_segmented.laz")

xyz = sf::st_coordinates(seeds)
seeds$Z = xyz[,3]

writeLAS(las, o)
terra::writeRaster(dtm, r)
#sf::st_write(sf::st_zm(seeds), s, append = FALSE)

# ==== EXPORT ALL TREES ====

trees_no_foliage = filter_poi(trees, foliage == FALSE)
plot(trees_no_foliage, color = "treeID", legend = TRUE, size = 2) |> add_dtm3d(dtm)
for (i in unique(trees_no_foliage$treeID))
{
  print(i)
  uid =  paste0(sample(c(letters, 0:9), 4, replace = TRUE), collapse = "")
  tree = filter_poi(trees_no_foliage, treeID == i)
  out = dirname(o)
  olas = paste0(out, "/ITS/tree_", i, "_", uid, ".las")
  oxyz = paste0(out, "/ITS/tree_", i, "_", uid, ".xyz")
  xyz = tree@data[, .(X,Y,Z)]
  writeLAS(tree, olas)
  data.table::fwrite(xyz, oxyz, sep = " ", col.names = FALSE)
}

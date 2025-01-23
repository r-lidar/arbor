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

# [MURRAY]: filte ris critical. Do no load the full point density but rather 10 to 20% only.
# I recommend trying to reach 15.000 20.000 pts/m². Load 20% of your point cloud. Here you shared
# a file decimated at 10% so I'm using filter = "" to load 100%
select = "tircn0"
filter = "-keep_random_fraction 0.2"
filter = ""
file = "~/Documents/Entreprise/clients/Forest Analysis Ltd/PRF/PRF/PRF025_15m_sor_10pct.laz"

# ====== READ POINT CLOUD =======

las = readTLS(file, select = select, filter = filter)
las

#las = clip_circle(las, 30,25,5)
#las = clip_circle(las, 70,-50,5)

# ====== REDUCE TO 1 POINTS/CM3 =======

# [MURRAY]: keeping one point per 1 cm voxel have actually barely no effect. It keeps almost all the points
# you can skip this step

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

# [MURRAY]. I'm working with 2 cm accurate point cloud. Yours is closer to 4 cm. I'm using 50. If
# you have more points or more noise increase the number
# Then plot to see if wood foliage looks roughly not too bad). Foliage should be blue
# Trunk and branches should be yellow and red. This is not the segmentation but this step should make sense
k_ani = 50

las = compute_anisotropy(las, k = k_ani)

if (display) plot(las, color = "anisotropy", legend = T, breaks = "quantile")

# ====== SEGMENT FOLIAGE/WOOD ======

las = segment_foliage(las, dtm, res = 0.08, max_gap = max_gap, k = k, min_passage = 5, th_anisotropy = 0.75)

if (display)
{
plot(las, color = "foliage", pal = c("#523b0a", "darkgreen"), size = 2) |> add_dtm3d(dtm)
plot(filter_poi(las, foliage == FALSE), pal = "#523b0a", size = 2) |> add_dtm3d(dtm)
x = plot(las, color = "foliage", pal = c("#523b0a", "darkgreen"), size = 2) |> add_dtm3d(dtm)
sk = filter_poi(las, skeleton == T)
plot(sk, add = x, pal = "red", size = 3)
}

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

# [MURRAY] Remove small tree that are likely badly segmented and clean the low understory
trees = clean_small_cluster(las, max_heigh = 2)

if (display)
{
plot(trees, color = "skeleton", legend = TRUE) |> add_dtm3d(dtm)
plot(trees, color = "treeID", legend = TRUE, size = 2) |> add_dtm3d(dtm)
plot(trees, color = "foliage", pal = c("chocolate4", "darkgreen")) |> add_dtm3d(dtm)
plot(filter_poi(trees, foliage == FALSE), color = "treeID", legend = TRUE) |> add_dtm3d(dtm) |> add_treetops3d(seeds, radius = 0.1)
}


# ==== CLIP BUFFER ======

# Remove the tree in the buffer

valid_trees = clip_buffer(trees, seeds, buffer = -2)
plot(valid_trees, color = "treeID", legend = TRUE, size = 2) |> add_dtm3d(dtm)|> add_dtm3d(dtm) |> add_treetops3d(seeds, radius = 0.1)

# ==== EXPORT TO FILE ======

o =  tools::file_path_sans_ext(file)
r = paste0(o, "_dtm.tif")
#s = paste0(o, "_seeds.shp")
o = paste0(o, "_segmented.laz")

xyz = sf::st_coordinates(seeds)
seeds$Z = xyz[,3]

writeLAS(las, o)
terra::writeRaster(dtm, r)
#sf::st_write(sf::st_zm(seeds), s, append = FALSE)

# ==== EXPORT ALL TREES ====

# To formats las and xyz

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

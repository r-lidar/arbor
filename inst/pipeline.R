library(lidR)
library(lidRtls)

z_factor = 0.8
k = 10
max_gap = 0.2
ncpu = 12
slice_seeds_at = c(0.75, 1.25) #; slice_seeds_at = c(0.8, 1)
path_finder_resolution = 0.1
voxel_sampling_resolution = 0.01
cut_above_ground = 0.5

select = "tircn0"
filter = "-keep_random_fraction 0.1" ; filter = ''
file = "~/Documents/Entreprise/clients/fsinvestor/SanDiego/FSTESTSCAN3UCSD_01_laz1_4_extract30m.laz"
file = "~/Documents/Entreprise/clients/fsinvestor/ST-X-ZamPlot/ZamPlot_part1_decimated.laz"
file = "~/Documents/Entreprise/clients/fsinvestor/ST-X-ZamPlot/ZamPlot_part2_decimated.laz"
file = "~/Documents/Entreprise/clients/fsinvestor/Rwanda/KwandahillsideHH_01_laz1_4_cliped_subsampled.laz"
file = "~/Téléchargements/PRF036_15m_sor.laz"
file = "~/Téléchargements/PRF025_15m_sor.laz"
#file = "~/Documents/Entreprise/clients/fsinvestor/Rwanda/Part3 body.laz"

# ====== READ POINT CLOUD =======

las = readTLS(file, select = select, filter = filter)
las

#las = clip_circle(las, 30,25,5)
#las = clip_circle(las, 70,-50,5)

# ====== REDUCE TO 1 POINTS/CM3 =======

algo = lowest_attribute_per_voxel(voxel_sampling_resolution, "Range")
las = decimate_points(las, algo)

las = lidR::knn_distance(las, k = 20)
plot(las, color = "distance", breaks = "quantile", legend = T)


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
mid = filter_poi(las, hag > cut_above_ground, hag < cut_above_ground + 3)
top = filter_poi(las, hag > cut_above_ground + 3)
las = filter_poi(las, hag > cut_above_ground)

plot(mid, color = "distance", breaks = "quantile", legend = T)
plot(top, color = "distance", breaks = "quantile", legend = T)
plot(las)

# ===== ANISOTROPY =======

las = compute_anisotropy(las, k = 50)
plot(las, color = "anisotropy", legend = T, breaks = "quantile")

# ====== SEGMENT FOLIAGE/WOOD ======

las = segment_foliage(las, dtm, max_gap = max_gap, k = k)

plot(las, color = "foliage", pal = c("#523b0a", "darkgreen"), size = 2) |> add_dtm3d(dtm)
plot(filter_poi(las, foliage == FALSE), pal = "#523b0a", size = 2) |> add_dtm3d(dtm)
x = plot(las, color = "foliage", pal = c("#523b0a", "darkgreen"), size = 2) |> add_dtm3d(dtm)
sk = filter_poi(las, skeleton == T)
plot(sk, add = x, pal = "red", size = 3)


#writeLAS(las, "~/Téléchargements/segmented.las")
#las = readLAS("~/Téléchargements/segmented.las", select = "xyz01234")

# ====== FIND TREE SEEDS =======

seeds = find_seeds(las, slice_seeds_at)

plot(las, color = "foliage", pal = c("chocolate4", "darkgreen")) |> add_treetops3d(seeds, radius = 0.1)
plot(filter_poi(las, hag < slice_seeds_at[2]), color = "foliage", pal = c("chocolate4", "darkgreen")) |> add_treetops3d(seeds, radius = 0.05)


# ====== SEGMENT TREES =======

las = segment_vegetation(las, seeds, max_gap = max_gap, k = k)

x = plot(las, color = "treeID") |> add_dtm3d(dtm) |> add_treetops3d(seeds, radius = 0.08)

trees = clean_small_cluster(las, max_heigh = 5)
plot(trees, color = "skeleton", legend = TRUE) |> add_dtm3d(dtm)
plot(trees, color = "treeID", legend = TRUE) |> add_dtm3d(dtm)
plot(trees, color = "foliage", pal = c("chocolate4", "darkgreen")) |> add_dtm3d(dtm) |> add_treetops3d(seeds, radius = 0.1)
plot(filter_poi(trees, foliage == FALSE), color = "treeID", legend = TRUE) |> add_dtm3d(dtm) |> add_treetops3d(seeds, radius = 0.1)

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
sf::st_write(sf::st_zm(seeds), s)

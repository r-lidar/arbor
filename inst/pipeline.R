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
file = "~/Documents/Entreprise/clients/fsinvestor/Rwanda/Forest site 1/Referencesite1_part2_subsample0.5.laz"
file = "~/Documents/Entreprise/clients/Forest Analysis Ltd/PRF/PRF/PRF025_15m_sor_10pct.laz" ; filter = ""
file = "~/Documents/Entreprise/clients/Forest Analysis Ltd/PRF/PRF/PRF193_15m_sor_10pct.laz" ; filter = ""
file = "~/Documents/Entreprise/clients/Forest Analysis Ltd/PRF/PRF/PRF200_15m_sor_10pct.laz" ; filter = ""

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

if (display) plot(las) |> add_dtm3d(dtm)

# ===== COMPUTE ANISOTROPY =======

#olas = sf::st_coordinates(las)
#'las = lidRtls:::smooth3d(las, 0.04)
#'las = lidR::knn_distance(las, k = 20)
#'f <- ecdf(las$distance)
#'las@data$anisotropy <- 1-f(las$distance)

las = compute_anisotropy(las, k = 50)

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

# ====== RETAIN ONLY MAIN TREES =======
#     (and clean the understory)

trees = clean_small_cluster(las, max_heigh = 5)
if (display) x = plot(trees, color = "treeID") |> add_dtm3d(dtm)

# ====== FIX SEGMENTATION ISSUES =======
# Segmentation is not always perfect especially complex environment.
# The seed dectection may have assign two seed to a single trees or an additionnal
# patch of wood may be assigned the id of a big tree because of a missing seed.
# For example in the example plot we have
if (display) plot(filter_poi(las, treeID %in% c(581,600)), color = "treeID")

trees = fix_split_trees(trees)
trees = fix_small_isolated_low_clusters(trees)

if (display) x = plot(trees, color = "treeID") |> add_dtm3d(dtm) |> add_treetops3d(seeds, radius = 0.08)

#' las$X = olas[,1]
#' las$Y = olas[,2]
#' las$Z = olas[,3]

if (display)
{
plot(trees, color = "skeleton", legend = TRUE) |> add_dtm3d(dtm)
plot(trees, color = "treeID", legend = TRUE, size = 2) |> add_dtm3d(dtm)
plot(trees, color = "foliage", pal = c("chocolate4", "darkgreen")) |> add_dtm3d(dtm)
plot(filter_poi(trees, foliage == FALSE), color = "treeID", legend = TRUE) |> add_dtm3d(dtm) |> add_treetops3d(seeds, radius = 0.1)
}

# ====== BUILD TREE EXTENSIONS =======
# (is a little slow, need improvement)

extensions = tree_extensions(trees, dtm, extra_height = 0.15)

if (display)
{
  x = plot(trees) |> add_dtm3d(dtm)
  plot(extensions, add = x)
}

trees = weld_extension(trees, extensions)

if (display)
{
  plot(trees, color = "treeID", legend = TRUE, size = 2) |> add_dtm3d(dtm)
  plot(trees, color = "foliage", pal = c("chocolate4", "darkgreen")) |> add_dtm3d(dtm)
}

# ==== CLIP BUFFER ======

valid_trees = clip_buffer(trees, seeds, -1)
plot(filter_poi(valid_trees, foliage == FALSE), color = "treeID", legend = TRUE) |> add_dtm3d(dtm) |> add_treetops3d(seeds, radius = 0.1)

# ==== VARIOUS EXPORTS ====

o =  tools::file_path_sans_ext(file)
r = paste0(o, "_dtm.tif")
s = paste0(o, "_seeds.shp")
o = paste0(o, "_segmented.laz")

xyz = sf::st_coordinates(seeds)
seeds$Z = xyz[,3]

writeLAS(las, o)
terra::writeRaster(dtm, r)
#sf::st_write(sf::st_zm(seeds), s, append = FALSE)

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

cmd = paste0("/home/jr/Logiciels/AdTree/Release/bin/AdTree " , out, "/ITS ", out ,"/QSM -radius 0.003 -alpha 0.8 -subtree 0.02")



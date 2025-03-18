library(lidR)
library(lidRtls)

set_lidr_threads(0)

foliage.colors = c("chocolate4", "darkgreen")

slice_seeds_at = c(0.75, 1)
cut_above_ground = 0.25

display = FALSE

select = ""


# [TIPS] The `filter` option is critical. Do not load the full point density, but rather only 20 to 30%.
# I recommend aiming for 20,000 pts/m². Loading more points makes everything slower but not
# necessarily more accurate. Of course 20.000 pts/m² means nothing in 3D and depends on the forest density
# but the idea is too drastically reduce the number of points.

file = "~/Documents/Entreprise/clients/fsinvestor/Rwanda/Kwandahillside/Kwandahillside.laz" ; filter = "-keep_random_fraction 0.3"
#file = "~/Documents/Entreprise/clients/fsinvestor/Rwanda/Kwandahillside/Kwandahillside_poisson.laz"

# ====== READ POINT CLOUD =======

# [TIPS] Do not use readLAS!

las = readTLS(file, select = select, filter = filter)

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
las@data$pointID = 1:npoints(las)

# ====== KEEP ABOVE DTM ======


bottom = filter_poi(las, hag <= cut_above_ground)
las = filter_poi(las, hag > cut_above_ground)

if (display) plot(las) |> add_dtm3d(dtm)


# ====== CLEAN BOTTOM NOISE  ======

# [TIPS] Removing noise only on the lowest layer to preserve foliage and high branches.
# Not recommended if the understory has a lot of foliage, a lot of sapling, a lot of mess!

if (FALSE)
{
  las = classify_noise(las, sor(m = 0.5))

  if (display)
  {
    plot(las, color = "Classification")
  }

  las = filter_poi(las, ! (Classification == LASNOISE & hag < 2))
}


# ===== COMPUTE ANISOTROPY =======

# [TIPS] Anisotropy is a critical stage. It is controlled mainly by the parameter k.
# k depends on the point density and sensor accuracy. More points? Increase k!
# Less accurate sensor? Increase k!
# k = 50 works well at 15,000–20,000 pts/m² with a sensor accuracy of 2–4 cm.
# Plot the output, then check if the wood foliage looks reasonable.
# Foliage should be blue, while trunks and branches should be yellow and red.
# This is not the actual segmentation, so it is not intended to be perfect,
# but this step should make sense and look roughly correct.


# Personal edge case. Not for public use
olas = sf::st_coordinates(las)
las = lidRtls:::smooth3d(las, 0.04)
las = lidR::knn_distance(las, k = 20)
f <- ecdf(las$distance)
las@data$anisotropy <- 1-f(las$distance)


if (display) plot(las, color = "anisotropy", legend = T, breaks = "quantile")

# ====== SEGMENT FOLIAGE/WOOD ======

# [TIPS] segment_foliage relies on a good anisotropy measurement. The method is described
# in the documentation

las = segment_foliage(las, dtm, res = 0.05, space_res = 0.1)

if (display)
{
  plot(las, color = "foliage", pal = foliage.colors, size = 2) |> add_dtm3d(dtm)
  plot(filter_poi(las, foliage == FALSE), pal = foliage.colors[1], size = 2) |> add_dtm3d(dtm)
  x = plot(las, color = "foliage", pal = foliage.colors, size = 2) |> add_dtm3d(dtm)
  sk = filter_poi(las, skeleton == T)
  plot(sk, pal = "red", size = 3)
}

# ====== FIND TREE SEEDS =======

# [TIPS] `find_seeds` is a critical step and probably the weak link in the pipeline.
# To segment individual trees, we need to identify seeds and assign an ID to each one.
# Then, a pathfinding algorithm determines the least-cost path from each point to a seed
# and assigns the ID of the nearest seed with the lowest cost.
#
# To find seeds, the algorithm takes a low slice, extracts points classified as wood,
# and identifies clusters. In scenes with dense understory vegetation, trunk obfuscation,
# noise, or lianas, this step can be particularly challenging.
#
# The most common issues are:
# 1. several seeds with different IDs for the same tree because the wood/foliage was not
#    perfectly segmented close to the ground or because of obfuscations of the trunks.
# 2. several trees have only one single seed because they are too close
#
# The most common issues are:
# 1. Several seeds with different IDs for the same tree because the wood/foliage was not
#    perfectly segmented near the ground or due to trunk obfuscation.
# 2. Multiple trees having only a single seed because they are too close together.
#
# For issue #1, there is a post-processing function available.
# For issue #2, this there is a reprocessing step

seeds = find_seeds(las, slice_seeds_at = slice_seeds_at)

if (display)
{
  col = pastel.colors(length(unique(seeds$treeID)))
  col = col[as.integer(as.factor(seeds$treeID))]
  plot(filter_poi(las,  hag > slice_seeds_at[1], hag < slice_seeds_at[2]), color = "foliage", pal = foliage.colors) |> add_treetops3d(seeds, radius = 0.05, color = col) |> add_dtm3d(dtm)
}

# ====== SEGMENT TREES =======

# [TIPS] `find_seeds` is the critical step here.

las = segment_vegetation(las, seeds, res = 0.05, max_gap = 0.2)

if (display)
{
  x = plot(las, color = "treeID") |> add_dtm3d(dtm)
}

# ====== RESTORE UNSMOOTHED COORDIANTES ======

las@data$X = olas[,1]
las@data$Y = olas[,2]
las@data$Z = olas[,3]

# ====== RETAIN ONLY MAIN TREES =======

# [TIPS] Low understory is unlikely to be properly segmented in complex contexts.
# In simpler contexts, this threshold can be set to a lower value.
# The goal is to retain the main trees and clean up the understory.

trees = clean_small_cluster(las, max_heigh = 3)

if (display) x = plot(trees, color = "treeID") |> add_dtm3d(dtm)

# ====== FIX SEGMENTATION ISSUES =======

# [TIPS] Segmentation is not always perfect, especially in complex environments.
# The seed detection may assign two seeds to a single tree, or an additional
# patch of wood may be assigned the ID of a large tree due to a missing seed.

trees = fix_small_isolated_low_clusters(trees)

# ==== CLIP BUFFER ======

if (FALSE)
{
  trees = clip_buffer(trees, seeds, -2)

  if (display)
  {
    plot(valid_trees, color = "treeID", legend = TRUE) |> add_dtm3d(dtm)
    plot(filter_poi(valid_trees, foliage == FALSE), color = "treeID", legend = TRUE) |> add_dtm3d(dtm)
  }
}

# ==== QSM ======

# QSM for a random tree

id = sample(unique(trees$treeID), 1)
tree = filter_poi(trees, treeID == id)
plot(tree, color = "foliage", pal = foliage.colors)

qsm <- qsm(tree, 1, 0.019, 0.002)
plot(tree, color = "foliage", pal = c("chocolate4", "darkgreen"), bg = "white", axis = T) |> add_qsm3d(qsm)

plot_qsm3d(qsm)
rgl::axes3d()

# Because it looks nice
x = plot_qsm3d(qsm, bottom_to_zero = FALSE)
plot(filter_poi(tree, foliage == TRUE), pal = "darkgreen", add = x, size = 2)
rgl::axes3d()


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

# ==== QSM ====

# Personal use only

cmd = paste0("/home/jr/Logiciels/AdTree/Release/bin/AdTree " , out, "/ITS ", out ,"/QSM -radius 0.003 -alpha 0.8 -subtree 0.02")
system(cmd)



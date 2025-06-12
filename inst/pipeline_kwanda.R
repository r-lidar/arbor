library(sf)
library(lidR)
library(lidRtls)

set_lidr_threads(0)

foliage.colors = c("chocolate4", "darkgreen")

slice_seeds_at = c(0.45, 0.65)
cut_above_ground = 0.25

display = FALSE

select = "i"


# [TIPS] The `filter` option is critical. Do not load the full point density, but rather only 20 to 30%.
# I recommend aiming for 20,000 pts/m². Loading more points makes everything slower but not
# necessarily more accurate. Of course 20.000 pts/m² means nothing in 3D and depends on the forest density
# but the idea is too drastically reduce the number of points.

file = "~/Documents/Entreprise/clients/fsinvestor/Rwanda/HageniaPlantation/Hagenia_2025-04-27_cliped_denoised.laz" ; filter = "-keep_random_fraction 0.25"

# ====== READ POINT CLOUD =======

# [TIPS] Do not use readLAS!

las = readTLS(file, select = select, filter = filter)

# ====== GROUND CLASSIFICATION ====

las = classify_ground(las, lidR::csf(rigidness = 3, class_threshold = 0.05, cloth_resolution = 0.05), last_returns = FALSE)
ground = filter_poi(las, Classification == LASGROUND)
ground = decimate_points(ground, lowest(0.25))
ground = classify_noise(ground, sor(k = 10, m = 2))
ground = remove_noise(ground)
ground$Classification = lidR::LASGROUND

# ====== DTM & HAG =======

dtm = rasterize_terrain(ground, 0.25, tin())
las = height_above_ground(las, algorithm = tin(), dtm = dtm)
las@data$pointID = 1:npoints(las)

# ====== KEEP ABOVE DTM ======

bottom = filter_poi(las, hag <= cut_above_ground)
las = filter_poi(las, hag > cut_above_ground)

if (display) plot(las) |> add_dtm3d(dtm)

chm = rasterize_canopy(las, 0.25)
las = merge_spatial(las, chm, "chm")
dz = las$Z - las$chm
las = filter_poi(las, Z - chm < -0.8)

# ====== CLEAN BOTTOM NOISE  ======

# [TIPS] Removing noise only on the lowest layer to preserve foliage and high branches.
# Not recommended if the understory has a lot of foliage, a lot of sapling, a lot of mess!

las = classify_noise(las, sor(m = 1))
plot(las, color = "Classification")

las = remove_noise(las)

#las = classify_noise(las, sor(k = 10, m = 0.25))

#if (display)
#{
#  plot(las, color = "Classification")
#}

#las = filter_poi(las, ! (Classification == LASNOISE & hag < slice_seeds_at[2]))



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
las = lidR::knn_distance(las, k = 10)
f <- ecdf(las$distance)
las@data$anisotropy <- 1-f(las$distance)

kkn = knn(las)
ann = las$anisotropy[kkn$nn.index]
ann = matrix(ann, ncol = ncol(kkn$nn.index),  nrow = nrow(kkn$nn.index))
ann = rowMeans(ann)
las@data$anisotropy = ann

if (display) plot(las, color = "anisotropy", legend = T, breaks = "quantile")

# ====== SEGMENT FOLIAGE/WOOD ======

# [TIPS] segment_foliage relies on a good anisotropy measurement. The method is described
# in the documentation

las = segment_foliage(las, dtm, res = 0.05, space_res = 0.1, min_passage = 2, th_anisotropy = 0.65)

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

las = segment_vegetation(las, seeds, res = 0.025, max_gap = 0.1)

if (display)
{
  x = plot(las, color = "treeID") |> add_dtm3d(dtm)
}

# ====== RESTORE UNSMOOTHED COORDIANTES ======

#as@data$X = olas[,1]
#las@data$Y = olas[,2]
#las@data$Z = olas[,3]

# ====== RETAIN ONLY MAIN TREES =======

# [TIPS] Low understory is unlikely to be properly segmented in complex contexts.
# In simpler contexts, this threshold can be set to a lower value.
# The goal is to retain the main trees and clean up the understory.

trees = clean_small_cluster(las, max_heigh = 2)

if (display) x = plot(trees[keep], color = "treeID") |> add_dtm3d(dtm)

# ====== FIX SEGMENTATION ISSUES =======

# [TIPS] Segmentation is not always perfect, especially in complex environments.
# The seed detection may assign two seeds to a single tree, or an additional
# patch of wood may be assigned the ID of a large tree due to a missing seed.

trees = fix_small_isolated_low_clusters(trees)

if (display)
{
  plot(trees, color = "treeID") |> add_dtm3d(dtm)
  plot(trees, color = "foliage", pal = foliage.colors) |> add_dtm3d(dtm)
}

seeds = find_seeds(trees, slice_seeds_at = slice_seeds_at-0.1)
trees = segment_vegetation(trees, seeds, res = 0.025, max_gap = 0.2)
trees = clean_small_cluster(las, max_heigh = 2)
trees = fix_small_isolated_low_clusters(trees)

if (display)
{
  plot(trees, color = "treeID") |> add_dtm3d(dtm)
  plot(trees, color = "foliage", pal = foliage.colors) |> add_dtm3d(dtm)
}


# ==== VARIOUS EXPORTS ====

las_export = TRUE
xyz_export = FALSE

o =  tools::file_path_sans_ext(file)
r = paste0(o, "_dtm.tif")
s = paste0(o, "_seeds.shp")
o = paste0(o, "_segmented.laz")
odir = dirname(o)
odir = paste0(odir, "/ITS/")
if (!dir.exists(odir)) dir.create(odir, recursive = TRUE)

writeLAS(las, o)
terra::writeRaster(dtm, r, overwrite = T)

trees_no_foliage = filter_poi(trees, foliage == FALSE)
plot(trees_no_foliage, color = "treeID", legend = TRUE, size = 2) |> add_dtm3d(dtm)
for (i in unique(trees_no_foliage$treeID))
{
  print(i)
  uid =  paste0(sample(c(letters, 0:9), 4, replace = TRUE), collapse = "")
  tree = filter_poi(trees_no_foliage, treeID == i)

  if (las_export)
  {
    olas = paste0(odir, "tree_", i,".las")
    writeLAS(tree, olas)
  }

  if (xyz_export)
  {
    oxyz = paste0(odir, "tree_", i,".xyz")
    xyz = tree@data[, .(X,Y,Z)]
    data.table::fwrite(xyz, oxyz, sep = " ", col.names = FALSE)
  }
}

# ==== QSM ====

# Personal use only

cmd = paste0("/home/jr/Logiciels/AdTree/Release/bin/AdTree " , out, "/ITS ", out ,"/QSM -radius 0.003 -alpha 0.8 -subtree 0.02")
system(cmd)



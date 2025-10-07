library(lidR)
library(lidRtls)

set_lidr_threads(0)

slice_seeds_at = c(0.5, 1)
cut_above_ground = 0.25

display = FALSE

# [TIPS] The `filter` option is critical. Do not load the full point density, but rather only 20 to 30%.
# I recommend aiming for 15,000 pts/m². Loading more points makes everything slower but not
# necessarily more accurate. Of course 20.000 pts/m² means nothing in 3D and depends on the forest density
# but the idea is too drastically reduce the number of points.

# California
file = "~/Documents/Entreprise/clients/fsinvestor/SanDiego/FSTESTSCAN3UCSD_01_laz1_4_extract30m.laz" ; filter = "-keep_random_fraction 0.3"

# Zambia
file = "~/Documents/Entreprise/clients/fsinvestor/Zambia/JasonHouse/ZamPlot_part1.laz" ; filter = "-keep_random_fraction 0.3"
file = "~/Documents/Entreprise/clients/fsinvestor/Zambia/JasonHouse/ZamPlot_part1_poisson.laz" ; filter = ""
file = "~/Documents/Entreprise/clients/fsinvestor/Zambia/JasonHouse/ZamPlot_part2.laz" ; filter = "-keep_random_fraction 0.3"
file = "~/Documents/Entreprise/clients/fsinvestor/Zambia/JasonHouse/ZamPlot_part3.laz" ; filter = "-keep_random_fraction 0.3"
file = "~/Documents/Entreprise/clients/fsinvestor/Zambia/JasonHouseTrees/Tree1-10_subsampled_rnd30.laz" ; filter = ""
file = "~/Documents/Entreprise/clients/fsinvestor/Zambia/JasonFarm/JasonFarm.laz" ; filter = ""
file = "~/Documents/Entreprise/clients/fsinvestor/Zambia/JasonFarm/JasonFarm_segmented.laz" ; filter = ""

# Indonesia
file = "~/Documents/Entreprise/clients/fsinvestor/Indonesia/Waykambawalk1RTK_01.las" ; filter = "-keep_random_fraction 0.3"
file = "~/Documents/Entreprise/clients/fsinvestor/Indonesia/Walk1Area2/Walk1area2slam_20x20plot.laz" ; filter = "-keep_random_fraction 0.3"
file = "~/Documents/Entreprise/clients/fsinvestor/Indonesia/MasiveButtroot/Waykananbindotree_01_isolated.laz" ; filter="-keep_random_fraction 0.6"

# Rwanda
file = "~/Documents/Entreprise/clients/fsinvestor/Rwanda/Kwandahillside/Kwandahillside.laz" ; filter = "-keep_random_fraction 0.2"
file = "~/Documents/Entreprise/clients/fsinvestor/Rwanda/Forest site 1/Referencesite1_part2_subsample0.5.laz" ; filter = ""
file = "~/Documents/Entreprise/clients/fsinvestor/Rwanda/Eucalyptus/Eucalyptuswalk2_01_subsampled_50.laz" ; filter = "-keep_random_fraction 0.333"

# Petawawa Research Forest
file = "~/Documents/Entreprise/clients/Forest Analysis Ltd/PRF/PRF/PRF025_15m_sor_10pct.laz" ; filter = ""
file = "~/Documents/Entreprise/clients/Forest Analysis Ltd/PRF/PRF/PRF193_15m_sor_10pct.laz" ; filter = ""
file = "~/Documents/Entreprise/clients/Forest Analysis Ltd/PRF/PRF/PRF200_15m_sor_10pct.laz" ; filter = ""
file = "~/Documents/Entreprise/clients/Forest Analysis Ltd/PRF/PRF/P0020_05_MLS_10m_buf10m_pj_z_range_30x30_test.las" ; filter = "-keep_random_fraction 0.3"

# Richard's data
file = "~/Documents/Usherbrooke/data/TN00/MLS-TN00-clip.laz" ; filter = "-keep_random_fraction 0.3"
file = "/home/jr/Documents/Usherbrooke/data/PRF002/MLS-PRF002-clip.laz" ; filter = "-keep_random_fraction 0.2"

# MRNF Oak plantations
file = "/home/jr/Documents/Entreprise/clients/MRNF-MLS/StA/test_plot1.las" ; filter = "-keep_random_fraction 0.8" ; cut_above_ground = 0.5


fdtm = paste0(tools::file_path_sans_ext(file), "_dtm.tif")

# ====== READ POINT CLOUD =======

# [TIPS] Do not use readLAS! use readTLS! It sorts the point cloud for cache
# efficiency

las = readTLS(file, select = "0", filter = filter)
las

plot(header(las))

# ===== GROUND CLASSIFICATION ======

las = classify_ground(las, lidR::csf(rigidness = 1, class_threshold = 0.05, cloth_resolution = 0.1), last_returns = FALSE)
ground = filter_poi(las, Classification == LASGROUND)
ground = decimate_points(ground, lowest(0.25))
ground = classify_noise(ground, sor(k = 10, m = 2))
ground = remove_noise(ground)
ground$Classification = lidR::LASGROUND

# ====== DTM ======

dtm = rasterize_terrain(ground, 0.5, tin())
terra::writeRaster(dtm, fdtm)

# ====== HEIGHT ABOVE GROUND =======

las = height_above_ground(las, algorithm = tin(), dtm = dtm)

# ====== KEEP ABOVE DTM ======

# [TIPS] We remove point close to the ground. It is impossible to segment
# anything close to the ground. This remove a lot of points, reduces computation time
# and clean the understory

bottom = filter_poi(las, hag <= cut_above_ground)
las = filter_poi(las, hag > cut_above_ground)

if (display) plot(las) |> add_dtm3d(dtm)


# ===== PRECOMPUTE DECIMATION =====

# [TIPS] In order to compute fast the point cloud must be decimated
# in some steps. The decimation technic is not straightforward and computationally
# demanding. Since the decimation is computed several times we can precompute it once
# and avoid recomputation by labeling retained points

las = barycentric_predecimation(las, 0.05)

if (display) plot(filter_poi(las, decimated == TRUE))

# ====== CLEAN BOTTOM NOISE  ======

# [TIPS] Removing noise only on the lowest layer to preserve foliage and high branches.
# Not recommended if the understory has a lot of foliage, a lot of sapling, a lot of mess!
# Might be useful is some specific case. Not recommanded.

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

las = compute_anisotropy(las, k = 50)

if (display) plot(las, color = "anisotropy", legend = T, breaks = "quantile")

# ====== SEGMENT FOLIAGE/WOOD ======

# [TIPS] segment_foliage relies on a good anisotropy measurement. The method is described
# in the documentation

las = segment_foliage(las, dtm, res = .05, min_passage = 5, max_gap = 1, k = 5)

if (display)
{
  plot(las, color = "foliage", pal = foliage.colors, size = 2) |> add_dtm3d(dtm)
  plot(filter_poi(las, foliage == FALSE), pal = foliage.colors[1], size = 2) |> add_dtm3d(dtm)
  x = plot(las, color = "foliage", pal = foliage.colors, size = 2) |> add_dtm3d(dtm)
  passage = filter_poi(las, passage > 0)
  passage@data$passage = log(passage$passage)
  plot(passage, color = "passage", legend = T)

  x = plot(filter_poi(las,  hag < 2.05), color = "foliage", pal = foliage.colors) |> add_dtm3d(dtm)
  plot(filter_poi(passage, hag< 2.05), add = x, legend = T, size = 4)

  x = plot(filter_poi(las,  (hag > 1 & hag < 1.05) | (hag > 2 &hag < 2.05), foliage == 0), color = "foliage", pal = foliage.colors) |> add_dtm3d(dtm)
  plot(filter_poi(passage, hag< 2.05), add = x, color = "passage", legend = T)
}

# ====== FIND TREE SEEDS =======

# [TIPS] `find_seeds` is a critical step and probably the weakest link in the pipeline.
# To segment individual trees, we need to identify seeds and assign an ID to each one.
# Then, a pathfinding algorithm determines the least-cost path from each point to a seed
# and assigns the ID of the seed with the lowest cost.
#
# To find seeds, the algorithm look at the previous paths taken during the foliage segmentation
# and aggregate then using connected component analysis.
#
# The most common issues are:
# 1. several seeds with different IDs for the same tree because the wood/foliage was not
#    perfectly segmented close to the ground or because of obfuscations of the trunks.
# 2. several trees have only one single seed because they are too close
#
# For issue #1, the new method should drastically reduced occurence of such case
# For issue #2, this there is a reprocessing step

somewood = filter_poi(las,  (hag > 1 & hag < 1.02) | (hag > 2 &hag < 2.02) | (hag > 3 & hag < 3.02), foliage == 0)
somewood = classify_noise(somewood, sor(k = 10, m = 0.5))
plot(somewood, color = "Classification")
somewood = remove_noise(somewood)
passages = filter_poi(las, passage > 0, hag < 3.1)
temp = rbind(somewood, passages)

temp$Z = temp$Z*0.5
temp = lidR::connected_components(temp, 0.1, 1, name = "treeID")
seeds = filter_poi(temp, passage > 0)
seeds$Z = seeds$Z/0.5

if (display)
{
  col = pastel.colors(length(unique(seeds$treeID)))
  col = col[as.integer(as.factor(seeds$treeID))]

  x = plot(somewood,  color = "foliage", pal = foliage.colors)  |> add_dtm3d(dtm)
  plot(seeds, color = "treeID", add = x, size = 2)
}

seeds = filter_poi(seeds,  (hag > 0.95 & hag < 1.05) | (hag > 1.95 & hag < 2.05), foliage == 0)
seeds = sf::st_as_sf(seeds)["treeID"]

if (display)
{
  col = pastel.colors(length(unique(seeds$treeID)))
  col = col[as.integer(as.factor(seeds$treeID))]

  x = plot(filter_poi(las, hag < 2), color = "foliage", pal = foliage.colors) |> add_treetops3d(seeds, radius = 0.08, color = col) |> add_dtm3d(dtm)
  plot(passage, add = x, pal = "gray", size = 2)
}

# ====== SEGMENT TREES =======

# [TIPS] `find_seeds` is the critical step here.

las = segment_vegetation(las, seeds, res = 0.05, max_gap = 0.6, k = 5)

if (display)
{
  x = plot(las, color = "treeID") |> add_dtm3d(dtm)
  plot(sk, add = x, pal = "red", size = 6)
  plot(sk)
}

# ====== RETAIN ONLY MAIN TREES =======

# [TIPS] Low understory is unlikely to be properly segmented in complex contexts.
# In simpler contexts, this threshold can be set to a lower value.
# The goal is to retain the main trees and clean up the understory.

trees = clean_small_cluster(las, max_heigh = 6)

if (display)
{
  plot(trees, color = "treeID", legend = TRUE) |> add_dtm3d(dtm)
  plot(trees, color = "foliage", pal = foliage.colors) |> add_dtm3d(dtm)
  plot(filter_poi(trees, foliage == FALSE), color = "treeID", legend = TRUE) |> add_dtm3d(dtm)
}

# ====== FIX SEGMENTATION ISSUES =======

# [TIPS] Segmentation is not always perfect, especially in complex environments.
# The seed detection may assign two seeds to a single tree, or an additional
# patch of wood may be assigned the ID of a large tree due to a missing seed.

trees = fix_small_isolated_low_clusters(trees)

if (display)
{
  plot(trees, color = "treeID", legend = TRUE) |> add_dtm3d(dtm)
  plot(trees, color = "foliage", pal = foliage.colors) |> add_dtm3d(dtm)
  plot(filter_poi(trees, foliage == FALSE), color = "treeID", legend = TRUE) |> add_dtm3d(dtm) |> add_treetops3d(seeds, radius = 0.1)
}


# ==== CLIP BUFFER ======

# Edge tree are necessarily bad by construction. We need to remove a buffer large enough
# to exclude trees connected to partial edge trees.

valid_trees = clip_buffer(trees, -4)

if (display)
{
  plot(valid_trees, color = "treeID", legend = TRUE) |> add_dtm3d(dtm)
  plot(filter_poi(valid_trees, foliage == FALSE), color = "treeID", legend = TRUE) |> add_dtm3d(dtm)
}

# ==== QSM ======

require(lidRqsm)

# QSM for a random tree

id = sample(unique(trees$treeID), 1)
tree = filter_poi(trees, treeID == id)

qsm <- lidRqsm::qsm_lidrqsm(tree)
xoff = qsm$startX[1]
yoff = qsm$startY[1]
qsm$startX = qsm$startX-xoff
qsm$startY = qsm$startY-yoff
qsm$endX = qsm$endX-xoff
qsm$endY = qsm$endY-yoff
x = plot(tree, color = "foliage", pal = c("chocolate4", "darkgreen"), bg = "white", axis = T, add = c(xoff, yoff))
lidRqsm::plot_qsm(qsm, add = c(0,0), color = "branch_order")

# Because it looks nice
if (FALSE)
{
  plot(filter_poi(tree, foliage == TRUE), pal = "darkgreen", add = c(xoff, yoff), size = 2)
  lidRqsm::plot_qsm(qsm, add = c(0, 0),  color = "branch_order")
  add_dtm3d( c(xoff, yoff), dtm)
}

# ==== VARIOUS EXPORTS ====

o =  tools::file_path_sans_ext(file)
s = paste0(o, "_seeds.shp")
t = paste0(o, "_trees.laz")
v = paste0(o, "_validtrees.laz")
o = paste0(o, "_segmented.laz")

writeLAS(las, o)
writeLAS(trees, t)
writeLAS(valid_trees, v)

trees_no_foliage = filter_poi(trees, foliage == FALSE)
plot(trees_no_foliage, color = "treeID", legend = TRUE, size = 2) |> add_dtm3d(dtm)
for (i in unique(trees_no_foliage$treeID))
{
  print(i)
  uid =  paste0(sample(c(letters, 0:9), 4, replace = TRUE), collapse = "")
  tree = filter_poi(trees_no_foliage, treeID == i)
  out = dirname(o)
  olas = paste0(out, "/ITS/tree_", i, "_", uid, ".las")
  writeLAS(tree, olas)
}



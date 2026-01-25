library(lidR)
library(arbor)

set_lidr_threads(0)

cut_above_ground = 0.25

display = FALSE

# The `filter` option is critical. **Do NOT** load the full point density. Use rather 20 to 30% max!
# I recommend targeting for 15.000 pts/m². Loading more points makes everything slower but not
# necessarily more accurate. Of course 15.000 pts/m² means nothing in 3D and depends on the forest density
# but the idea is too drastically reduce the number of points from the original point cloud. Working
# full density is worthless. It blows up RAM consumption, blows up computation time and brings no value.

# California
file = "~/Documents/r-lidar/clients/fsinvestor/SanDiego/FSTESTSCAN3UCSD_01_laz1_4_extract30m.laz" ; filter = "-keep_random_fraction 0.3"

# Zambia
file = "~/Documents/r-lidar/clients/fsinvestor/Zambia/JasonHouse/ZamPlot_part1.laz" ; filter = "-keep_random_fraction 0.3"
file = "~/Documents/r-lidar/clients/fsinvestor/Zambia/JasonHouse/ZamPlot_part2.laz" ; filter = "-keep_random_fraction 0.15"
file = "~/Documents/r-lidar/clients/fsinvestor/Zambia/JasonHouse/ZamPlot_part3.laz" ; filter = "-keep_random_fraction 0.3"
file = "~/Documents/r-lidar/clients/fsinvestor/Zambia/JasonHouseTrees/Tree1-10_subsampled_rnd30.laz" ; filter = ""
file = "~/Documents/r-lidar/clients/fsinvestor/Zambia/JasonFarm/JasonFarm.laz" ; filter = ""
file = "~/Documents/r-lidar/clients/fsinvestor/Zambia/JasonFarm/JasonFarm_segmented.laz" ; filter = ""

# Indonesia
file = "~/Documents/r-lidar/clients/fsinvestor/Indonesia/Waykambawalk1RTK_01.las" ; filter = "-keep_random_fraction 0.3"
file = "~/Documents/r-lidar/clients/fsinvestor/Indonesia/Walk1Area2/Walk1area2slam_20x20plot.laz" ; filter = "-keep_random_fraction 0.3"
file = "~/Documents/r-lidar/clients/fsinvestor/Indonesia/MasiveButtroot/Waykananbindotree_01_isolated.laz" ; filter="-keep_random_fraction 0.6"

# Rwanda
file = "~/Documents/r-lidar/clients/fsinvestor/Rwanda/Kwandahillside/Kwandahillside.laz" ; filter = "-keep_random_fraction 0.2"
file = "~/Documents/r-lidar/clients/fsinvestor/Rwanda/Forest site 1/Referencesite1_part2_subsample0.5.laz" ; filter = ""
file = "~/Documents/r-lidar/clients/fsinvestor/Rwanda/Eucalyptus/Eucalyptuswalk2_01_subsampled_50.laz" ; filter = "-keep_random_fraction 0.333"
file = "~/Documents/r-lidar/clients/fsinvestor/Rwanda/Eucalyptus/Stoneucsd6_01_subsampled_clio_40x40.laz" ; filter = ""


# Petawawa Research Forest
file = "~/Documents/r-lidar/clients/Forest Analysis Ltd/PRF/PRF/PRF025_15m_sor_10pct.laz" ; filter = ""
file = "~/Documents/r-lidar/clients/Forest Analysis Ltd/PRF/PRF/PRF193_15m_sor_10pct.laz" ; filter = ""
file = "~/Documents/r-lidar/clients/Forest Analysis Ltd/PRF/PRF/PRF200_15m_sor_10pct.laz" ; filter = ""
file = "~/Documents/r-lidar/clients/Forest Analysis Ltd/PRF/PRF/P0020_05_MLS_10m_buf10m_pj_z_range_30x30_test.las" ; filter = "-keep_random_fraction 0.3"

# Richard's data
file = "/home/jr/Documents/Usherbrooke/Registration/data/TN00/MLS-TN00-clip.laz" ; filter = "-keep_random_fraction 0.25"
file = "/home/jr/Documents/Usherbrooke/Registration/data/PRF002/MLS-PRF002-clip.laz" ; filter = "-keep_random_fraction 0.1"

# MRNF Oak plantations
file = "/home/jr/Documents/r-lidar/clients/MRNF-MLS/las/Tutoriel/sta_plot1_20x20.laz" ; filter = "-keep_random_fraction 0.8" ; cut_above_ground = 0.5

# Bastien's data
file = "~/Téléchargements/GJ-019_plot_15m_prep.las" ; filter = "-keep_random_fraction 0.25"
file = "~/Téléchargements/P2_clean.laz" ; filter = "-keep_random_fraction 0.2"
file = "~/Téléchargements/P1_clean_subset.laz" ; filter = "-keep_random_fraction 0.4"
file = "~/Téléchargements/P0024_07_segmented_subset.laz" ; filter = ""
file = "~/Téléchargements/P1_decimated.laz" ; filter = "-keep_random_fraction 0.6"

# French Guyana'
file = "/home/jr/Documents/r-lidar/clients/IRF/TropiscatNE/TropiscatNE_01_laz1_4_decimated_40x40.laz" ; filter = "-keep_random_fraction 0.4" ; cut_above_ground = 0.5


# UAV
file = "/home/jr/Téléchargements/Forestertestzone.las" ; filter = "" ; cut_above_ground = 0.25

# ===== PROCESSING PARAMETERS =====

params = default_arbor_parameters

# ====== READ POINT CLOUD =======

# Do not use readLAS use readTLS! It sorts the point cloud for L1 cache efficiency

las <- readTLS(file, select = "xyz", filter = filter)
las <- hybrid_homogeneization(las)

plot(header(las))

#las = clip_circle(las, mean(las$X), mean(las$Y), 15)

# Print to see your density. Target is between 10.000 and 20.000
print(las)
plot(header(las))

# ===== GROUND CLASSIFICATION ======

las    <- lidR::classify_ground(las, lidR::csf(rigidness = 1, class_threshold = 0.05, cloth_resolution = 0.1), last_returns = FALSE)
ground <- lidR::filter_poi(las, Classification == LASGROUND)
ground <- lidR::decimate_points(ground, lowest(0.25))
ground <- lidR::classify_noise(ground, sor(k = 10, m = 2))
ground <- lidR::remove_noise(ground)
ground$Classification <- lidR::LASGROUND
gc()

# ====== DTM & HAG ======

dtm <- lidR::rasterize_terrain(ground, 0.5, lidR::tin())
las <- lidR::height_above_ground(las, algorithm = lidR::tin(), dtm = dtm)

if (display) plot_dtm3d(dtm)

# ====== KEEP ABOVE DTM ======

# We remove points close to the ground. It is impossible to segment
# anything close to the ground. This remove a lot of points, reduces computation time
# and clean the understory. This is controlled by cut_above_ground. The value must
# be chosen depending on the level of understory complexity close to the ground. The idea
# is to remove most of the very low vegetation. 25 cm might be good. Some plot require 50 cm.

las <- lidR::filter_poi(las, hag > cut_above_ground)

if (display) plot(las) |> add_dtm3d(dtm)

gc()

# ===== PRECOMPUTE DECIMATION =====

# In order to compute fast the point cloud must be temporarily decimated during the computation
# in some steps. The decimation technic is not straightforward (no voxel decimation) and a little
# computationally demanding. Since the decimation is computed several times we can precompute it once
# and avoid recomputation by labeling retained points. Default 5 cm decimation is good.

las <- barycentric_predecimation(las, params)

if (display) plot(filter_poi(las, decimated == TRUE))

# ====== CLEAN BOTTOM NOISE  ======

# Removing noise only on the lowest layer to preserve foliage and high branches.
# Not recommended if the understory has a lot of foliage, a lot of sapling, a lot of mess!
# Might be useful is some specific case. Not recommanded to apply.

if (FALSE)
{
  las <- classify_noise(las, sor(m = 0.5))
  if (display) plot(las, color = "Classification")
  las <- filter_poi(las, !(Classification == LASNOISE & hag < 2))
}


# ===== COMPUTE ANISOTROPY =======

# Anisotropy is a critical stage. It is controlled mainly by the parameter k for knn search.
# Choice of k depends on the point density and sensor accuracy. More points? Increase k!
# Less accurate sensor? Increase k!
# k = 50 works well at 15.000–20.000 pts/m² with a sensor accuracy of 2–4 cm.
# Plot the output, then check if the wood foliage looks reasonable.
# Foliage should be blue, while trunks and branches should be yellow and red.
# This is not the actual segmentation, so it is not intended to be perfect,
# but this step should make sense and look roughly correct.

las <- compute_anisotropy(las, params)

if (display) plot_anisotropy(las)

# ====== SEGMENT FOLIAGE/WOOD ======

# segment_foliage relies, at least partially, on a good anisotropy measurement.
# If the previous step is bad this step will be bad too.

las <- segment_foliage(las, dtm, params)

if (display)
{
  plot_semantic(las, dtm) # Wood/foliage
  plot(filter_poi(las, foliage == FALSE), pal = foliage.colors[1], size = 2) |> add_dtm3d(dtm)   # Wood only
  plot_passage(las, dtm, th = 3)   # Pathfinder passages

  passage <- lidR::filter_poi(las, passage > 3)
  x <- plot_semantic(las, dtm)
  plot(passage, add = x, legend = T, size = 4)

  # Pathfinder passages + scene < 2m
  passage <- lidR::filter_poi(las, passage <15)
  x <- plot(filter_poi(las,  hag < 3), color = "foliage", pal = foliage.colors) |> add_dtm3d(dtm)
  plot_passage(passage, add = x, size = 3, th = 15)
}

# ====== FIND TREE SEEDS =======

# Finding seeds is a critical step and probably the weakest link in the pipeline.
# To segment individual trees, we need to identify seeds and assign an ID to each one.
# Then, a pathfinder determines the least-cost path from each point to a seed
# and assigns the ID of the seed with the lowest cost.
#
# To find seeds, the algorithm look at the previous paths taken during the foliage segmentation
# and aggregate them using connected component analysis.

params$path_finder$max_gap = 1

seeds <- find_seeds(las, params)

if (display)
{
  x <- plot(lidR::filter_poi(las,  hag < 2), color = "foliage", pal = foliage.colors, size =2) |> add_dtm3d(dtm)
  plot(seeds, color = "treeID", add = x, size = 8)
}

# ====== SEGMENT TREES =======

# Finding seed is THE critical step here.

las <- segment_vegetation(las, seeds, params)

if (display)
{
  x = plot_instance(las, dtm)
  plot(seeds, color = "treeID", add = x, size = 8)
  plot_semantic_instance(las, dtm)

}

# ====== RETAIN ONLY MAIN TREES =======

# Low understory is unlikely to be properly segmented in complex contexts.
# In simpler contexts, this threshold can be set to a lower value.
# The goal is to retain the main trees and clean up the understory. It also
# remove blob of points with no ID

trees <- remove_small_trees(las, max_height = 4)

if (display)
{
  plot_instance(trees, dtm)
  plot_semantic(trees, dtm)
  plot_semantic_instance(trees, dtm)
  plot(lidR::filter_poi(trees, foliage == FALSE), color = "treeID", legend = TRUE) |> add_dtm3d(dtm)
}

# ====== FIX SEGMENTATION ISSUES =======

# Segmentation is not always perfect, especially in complex environments.
# The seed detection may assign two seeds to a single tree, or an additional
# patch of wood may be assigned the ID of a large tree due to a missing seed.

trees <- fix_small_isolated_low_clusters(trees)

if (display)
{
  plot_semantic_instance(trees, dtm)
  plot_semantic(trees, dtm)
}

# ==== CLIP BUFFER ======

# Edge trees are necessarily bad by construction. We need to remove a buffer large enough
# to exclude trees connected to partially visible edge trees.

valid_trees <- clip_buffer(trees, -5)

if (display)
{
  plot(valid_trees, color = "treeID", legend = TRUE) |> add_dtm3d(dtm)
}

# ==== VARIOUS EXPORTS ====

o <- tools::file_path_sans_ext(file)
t <- paste0(o, "_trees.laz")
v <- paste0(o, "_validtrees.laz")
o <- paste0(o, "_segmented.laz")
r <- paste0(o, "_dtm.tif")

las = colorize_trees(las)
trees = colorize_trees(trees)
valid_trees = colorize_trees(valid_trees)


writeLAS(las, o)
writeLAS(trees, t)
writeLAS(valid_trees, v)
terra::writeRaster(dtm, r)

plot(valid_trees, color = "treeID", legend = TRUE, size = 2) |> add_dtm3d(dtm)
for (i in unique(valid_trees$treeID))
{
  print(i)

  tree <- filter_poi(valid_trees, treeID == i)
  olas <- paste0(dirname(o), "/ITS/tree_", i, ".las")
  writeLAS(tree, olas)
}
# ==== QSM ======

# QSM for a random tree
id   <- sample(unique(trees$treeID), 1)
tree <- lidR::filter_poi(trees, treeID == id)
plot_semantic(tree)
qsm  <- qsm(tree, step = 0.1, cl_dist = 0.2)
x <- plot_semantic(tree)
plot_qsm(qsm, color = "branch_order", add = x + c(-5, 0), skeleton = F)

# ==== QSM batch =====

qsm_batch( paste0(dirname(o), "/ITS/"), paste0(dirname(o), "/qsm/"))

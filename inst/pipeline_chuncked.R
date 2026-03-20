library(lidR)
library(arbor)

file = "/home/jr/Documents/r-lidar/clients/geotree/TLS/wytham/wytham_benchmark.las"
file = "/home/jr/Documents/r-lidar/clients/geotree/TLS/robsoncreek/robsoncreek_benchmark.las"
file = "/home/jr/Documents/r-lidar/clients/geotree/TLS/ofental/ofental_benchmark.las"
file = "/home/jr/Documents/r-lidar/clients/geotree/TLS/litchfield/litchfield_benchmark.las"

filter =  ""
cut_above_ground = 0.15

o <- tools::file_path_sans_ext(file)
dir <- paste0(o, "_output")
o = basename(o)
dir.create(dir, F)
tmp = paste0(dir, "/", "tmp")
dir.create(tmp, F)

ctg = readTLScatalog(file, select = "xyz0", filter = filter)
opt_chunk_size(ctg) = 50
opt_chunk_buffer(ctg) = 10
plot(ctg, chunk = T)
chunks = lidR:::engine_chunks(ctg)

params = default_arbor_parameters
params$path_finder$max_gap = 1
params$path_finder$k_neighborhood_connectivity = 20

display = FALSE

maxID = 0
for (i in seq_along(chunks))
{
  cat("==============\n")
  cat("Chunk", i, "of", length(chunks), "\n")
  cat("==============\n")

  f = chunks[[i]]
  las <- readLAS(f)
  las <- hybrid_homogeneization(las)

  # ===== GROUND CLASSIFICATION ======

  las    <- lidR::classify_ground(las, lidR::csf(rigidness = 1, class_threshold = 0.05, cloth_resolution = 0.1), last_returns = FALSE)
  ground <- lidR::filter_poi(las, Classification == LASGROUND)
  ground <- lidR::decimate_points(ground, lidR::lowest(0.25))
  ground <- lidR::classify_noise(ground, lidR::sor(k = 10, m = 2))
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

  # ===== COMPUTE ANISOTROPY =======

  # Anisotropy is a critical stage. It is controlled mainly by the parameter k for knn search.
  # Choice of k depends on the point density and sensor accuracy. More points? Increase k!
  # Less accurate sensor? Increase k!
  # k = 50 works well at 15.000–20.000 pts/m² with a sensor accuracy of 2–4 cm.
  # Plot the output, then check if the wood foliage looks reasonable.
  # Foliage should be blue, while trunks and branches should be yellow and red.
  # This is not the actual segmentation, so it is not intended to be perfect,
  # but this step should make sense and look roughly correct.

  las <- wood_likelihood(las, params)

  if (display) plot_likelihood(las)

  # ====== SEGMENT FOLIAGE/WOOD ======

  # segment_foliage relies, at least partially, on a good anisotropy measurement.
  # If the previous step is bad this step will be bad too.

  las <- segment_semantic(las, dtm, params)

  if (display)
  {
    plot_semantic(las, dtm) # Wood/foliage
  }

  # ====== FIND TREE SEEDS =======

  # Finding seeds is a critical step and probably the weakest link in the pipeline.
  # To segment individual trees, we need to identify seeds and assign an ID to each one.
  # Then, a pathfinder determines the least-cost path from each point to a seed
  # and assigns the ID of the seed with the lowest cost.
  #
  # To find seeds, the algorithm look at the previous paths taken during the foliage segmentation
  # and aggregate them using connected component analysis.

  seeds <- find_seeds(las, params)

  if (display)
  {
    x <- plot(lidR::filter_poi(las,  hag < 2), color = "foliage", pal = foliage.colors, size =2) |> add_dtm3d(dtm)
    plot(seeds, color = "treeID", add = x, size = 8)
  }

  # ====== SEGMENT TREES =======

  # Finding seed is THE critical step here.

  las <- segment_instance(las, seeds, params)

  if (display)
  {
    x = plot_instance(las, dtm)
    plot(seeds, color = "treeID", add = x, size = 8)
    plot_semantic_instance(las, dtm)
  }

  bb = st_bbox(f)
  bb = sf::st_as_sfc(bb)
  tm = seeds@data[, .(x = mean(X), y = mean(Y)), by = treeID]
  tm = sf::st_as_sf(tm, coords = c("x", "y"))
  seeds_roi = sf::st_contains(bb, tm)
  seeds_roi = tm[seeds_roi[[1]], ]
  seeds_roi_ids = unique(seeds_roi$treeID)
  seeds_roi = filter_poi(seeds, treeID %in% seeds_roi_ids)

  las = filter_poi(las, treeID %in% seeds_roi_ids)
  las@data$treeID = las@data$treeID + maxID
  maxID = max(las@data$treeID)

  dtm = terra::crop(dtm, terra::vect(bb))


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

  trees <- remove_small_trees(las, max_height = 2)

  if (display)
  {
    plot_instance(trees, dtm)
    plot_semantic(trees, dtm)
    plot_semantic_instance(trees, dtm)
    plot(lidR::filter_poi(trees, foliage == FALSE), color = "treeID", legend = TRUE) |> add_dtm3d(dtm)
  }

  # ==== VARIOUS EXPORTS ====

  t <- paste0(tmp, "/", o, "_trees_part", i, ".laz")
  s <- paste0(tmp, "/", o, "_segmented_part", i, ".laz")
  r <- paste0(tmp, "/", o, "_dtm_part", i, ".tif")

  las = colorize_trees(las, darken_foliage = F)
  trees = colorize_trees(trees, darken_foliage = F)

  writeLAS(las, s)
  writeLAS(trees, t)
  terra::writeRaster(dtm, r)
}

gc()

files = list.files(tmp, pattern = "dtm", full.names = T)
dtm = terra::vrt(files)
terra::writeRaster(dtm, paste0(dir, "/", o, "_dtm.tif"))

ctg = readLAScatalog(tmp, pattern = "segmented", select = "xyzic0")
las = readLAS(ctg)
writeLAS(las, paste0(dir, "/", o, "_full.laz"))
las = clip_buffer(las)
writeLAS(las, paste0(dir, "/", o, "_nobuffer.laz"))
rm(las)
gc()

ctg = readLAScatalog(tmp, pattern = "trees", select = "xyzic0")
las = readLAS(ctg)
las = clip_buffer(las)
writeLAS(las, paste0(dir, "/", o, "_trees.laz"))
rm(las)
gc()

library(lidR)
library(arbor)

file = "/home/jr/Documents/r-lidar/clients/geotree/TLS/wytham/wytham_benchmark.las"
file = "/home/jr/Documents/r-lidar/clients/geotree/TLS/robsoncreek/robsoncreek_benchmark.las"
file = "/home/jr/Documents/r-lidar/clients/geotree/TLS/ofental/ofental_benchmark.las"
file = "/home/jr/Documents/r-lidar/clients/geotree/TLS/litchfield/litchfield_benchmark.las"
file = "/home/jr/Documents/r-lidar/clients/geotree/TLS/bicuar/bicuar.P02.las"

display = FALSE
random_fraction  = 0.25
cut_above_ground = 0.15
buffer_size = 10

par = arbor_parameters_default
par$global$cut_above_ground = cut_above_ground

o <- tools::file_path_sans_ext(file)
dir <- paste0(o, "_output")
o = basename(o)
tmp = paste0(dir, "/", "tmp")
dir.create(dir, showWarnings = FALSE)
dir.create(tmp, showWarnings = FALSE)

ctg = readTLScatalog(file, select = "xyz0", filter = paste("-keep_random_fraction", random_fraction))
opt_chunk_size(ctg) = 50
opt_chunk_buffer(ctg) = buffer_size
plot(ctg, chunk = T)
chunks = lidR::engine_chunks(ctg)

maxID = 0
for (i in seq_along(chunks))
{
  cat("==============\n")
  cat("Chunk", i, "of", length(chunks), "\n")
  cat("==============\n")

  # Regular pipeline
  chk <- chunks[[i]]
  las <- readLAS(chk)
  las <- hybrid_homogeneization(las)
  las <- segment_ground(las, par)
  las <- wood_likelihood(las, par)
  las <- segment_semantic(las, par)
  see <- find_seeds(las, par)
  las <- segment_instance(las, see, par)
  las <- colorize_trees(las, FALSE)

  # Custom tree remove
  bb <- st_bbox(chk)
  bb <- sf::st_as_sfc(bb)
  tm <- see@data[, .(x = mean(X), y = mean(Y)), by = treeID]
  tm <- sf::st_as_sf(tm, coords = c("x", "y"))
  seeds_roi <- sf::st_contains(bb, tm)
  seeds_roi <- tm[seeds_roi[[1]], ]
  seeds_roi_ids <- unique(seeds_roi$treeID)
  seeds_roi <- filter_poi(seeds, treeID %in% seeds_roi_ids)

  las <- filter_poi(las, treeID %in% seeds_roi_ids)
  las$treeID <- las$treeID + maxID
  maxID <- max(las@data$treeID)

  dtm <- rasterize_terrain(las, 0.1)
  dtm <- terra::crop(dtm, terra::vect(bb))

  # Export
  s <- paste0(tmp, "/", o, "_segmented_part", i, ".laz")
  r <- paste0(tmp, "/", o, "_dtm_part", i, ".tif")
  writeLAS(las, s)
  terra::writeRaster(dtm, r)

  # Clean up
  rm(las)
  rm(see)
  rm(dtm)
  gc()
}

gc()

# Build the full DTM
files = list.files(tmp, pattern = "dtm", full.names = T)
dtm = terra::vrt(files)
terra::writeRaster(dtm, paste0(dir, "/", o, "_dtm.tif"))

# Build the full point cloud without buffer
ctg = readLAScatalog(tmp, pattern = "segmented", select = "xyzic0")
las = readLAS(ctg)
las = flag_buffer(las, buffer = -buffer_size)
las = filter_poi(las, treeID > 0)
writeLAS(las, paste0(dir, "/", o, "_full.laz"))

q("no")

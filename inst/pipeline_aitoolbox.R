split_rectangle_into_squares <- function(polygon, n)
{
  # Ensure it's projected (not in lon/lat) for square distance units
  if (sf::st_is_longlat(polygon)) {
    stop("Reproject the polygon to a projected CRS (e.g., UTM) before using this function.")
  }

  # Bounding box
  bbox <- st_bbox(polygon)
  width <- bbox$xmax - bbox$xmin
  height <- bbox$ymax - bbox$ymin

  # Total area
  total_area <- as.numeric(st_area(polygon))
  target_square_area <- total_area / n
  target_square_size <- sqrt(target_square_area)

  # Number of columns and rows
  ncol <- ceiling(width / target_square_size)
  nrow <- ceiling(height / target_square_size)

  # Build grid
  grid <- sf::st_make_grid(polygon, n = c(ncol, nrow), square = TRUE)

  # Intersect with original polygon to clean up edges
  grid <- sf::st_intersection(grid, polygon)

  return(grid)
}


library(lidR)
library(lidRtls)
library(sf)

set_lidr_threads(0)

foliage.colors = c("chocolate4", "darkgreen")

display = FALSE

file = "~/Documents/Entreprise/clients/fsinvestor/Zambia/JasonFarm/JasonFarm_segmented.laz" ; filter = ""

fdtm = paste0(tools::file_path_sans_ext(file), "_dtm.tif")

ctg = readLAScatalog(file, select = "0")
bbox = st_as_sf(ctg)
bbbox = st_as_sfc(st_buffer(bbox, -8, endCapStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 10))

chunks = split_rectangle_into_squares(bbox, 4)

plot(ctg)
plot(chunks, add = T)
plot(bbbox, add = T, lty = 3, border = "red")


for (i in seq_along(chunks))
{
  cat("===============\nChunk", i, "of", length(chunks), "\n===============\n\n")
  chunk = st_as_sf(chunks[i, ])
  bchunk = st_buffer(chunk, 8, endCapStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 10)

  if (display)
  {
    plot(ctg)
    plot(chunks, add = T)
    plot(chunk, add = T, col = 'lightblue')
    plot(bchunk, add = T, lty = 3)
    plot(bbbox, add = T, lty = 3, border = "red")
  }

  las = clip_roi(ctg, bchunk)

  if (display)
  {
    plot(las)
  }

  if (!file.exists(fdtm)) {

    # ---- GROUND CLASSIFICATION ----

    las = classify_ground(las, lidR::csf(rigidness = 1, class_threshold = 0.05, cloth_resolution = 0.1), last_returns = FALSE)
    ground = filter_poi(las, Classification == LASGROUND)
    ground = decimate_points(ground, lowest(0.25))
    ground = classify_noise(ground, sor(k = 10, m = 2))
    ground = remove_noise(ground)
    ground$Classification = lidR::LASGROUND

    # ---- DTM ---

    dtm = rasterize_terrain(ground, 0.5, tin())
    terra::writeRaster(dtm, fdtm)
  } else {
    dtm = terra::rast(fdtm)
  }

  dtm = terra::crop(dtm, bchunk)

  # ====== HEIGHT ABOVE GROUND =======

  las = height_above_ground(las, algorithm = dtm)

  if (display) plot(las) |> add_dtm3d(dtm)

  # ===== COMPUTE ANISOTROPY =======

  las = compute_anisotropy(las, k = 50)

  if (display) plot(las, color = "anisotropy", legend = T, breaks = "quantile")

  # ====== SEGMENT FOLIAGE/WOOD ======

  # [TIPS] segment_foliage relies on a good anisotropy measurement. The method is described
  # in the documentation

  las = segment_foliage(las, dtm, res = .05, min_passage = 5, max_gap = 1)

  las = add_lasattribute(las, las$itc, name = "treeID", desc = "Individual Tree ID")
  las = add_lasattribute(las, name = "hag", desc = "Height Above Ground")
  las = remove_lasattribute(las, "treefilter")
  las = remove_lasattribute(las, "itc")

  if (display)
  {
    plot(las, color = "foliage", pal = foliage.colors, size = 2) |> add_dtm3d(dtm)
    plot(filter_poi(las, foliage == FALSE), pal = foliage.colors[1], size = 2) |> add_dtm3d(dtm)
    x = plot(las, color = "foliage", pal = foliage.colors, size = 2) |> add_dtm3d(dtm)
    sk = filter_poi(las, skeleton == T)
    plot(sk, pal = "red", size = 3)
  }

  locate_seed = function(x,y,z)
  {
    i = order(z)
    x = x[i][1:min(length(x),100)]
    y = y[i][1:min(length(x),100)]
    list(x = mean(x), y = mean(y))
  }

  seeds = las@data[, locate_seed(X,Y,Z), by = treeID]
  seeds = st_as_sf(seeds, coords = c("x", "y"), crs = st_crs(las))

  valid_seeds = sf::st_filter(seeds, chunk)
  valid_seeds = sf::st_filter(valid_seeds, bbbox)

  plot(ctg)
  plot(chunks, add = T)
  plot(chunk, add = T, col = 'lightblue')
  plot(bchunk, add = T, lty = 3)
  plot(seeds, add = T, cex = 0.3, col = "red", pch = 19)
  plot(valid_seeds, add = T, cex = 0.3, col = "darkgreen", pch = 19)

  valid_trees = las[las$treeID %in% valid_seeds$treeID]
  valid_trees = clean_small_cluster(valid_trees, max_heigh = 4)
  idx = which(valid_trees@data$stemcls == 2)
  valid_trees@data$foliage[idx] = 0

  # ==== VARIOUS EXPORTS ====

  o =  tools::file_path_sans_ext(file)

  if (display) plot(valid_trees, color = "foliage", pal = foliage.colors,  legend = TRUE, size = 2) |> add_dtm3d(dtm)

  for (j in unique(valid_trees$treeID))
  {
    #uid =  paste0(sample(c(letters, 0:9), 4, replace = TRUE), collapse = "")
    tree = filter_poi(valid_trees, treeID == j)
    out = dirname(o)
    olas = paste0(out, "/its/tree_", j, ".las")
    writeLAS(tree, olas)
  }
}


cmd_blender <- function(args) {

  # --- Segment Usage ---
  usage_blender <- function() {
    cat("
Usage:
  arbor blender <input.las> [options]

Mandatory:
  <input.las>             Input LAS/LAZ file

Options:
  -cut <m>                Height threshold above ground [default: 0.25]
  -height <m>             Remove trees smaller than this [default: 2]
  -buffer <m>             Buffer removed from borders [default: 0]
  -fraction <0-1>         Keep random fraction [default: 0.25]

Other:
  -h, --help               Show this help
")
    quit(save = "no", status = 0)
  }

  if (length(args) == 0 || has_flag(args, "-h") || has_flag(args, "--help")) {
    usage_segment()
  }

  # --- Parse Inputs ---
  input <- args[1]
  if (startsWith(input, "-")) {
    fail("First argument must be the input file.")
  }
  if (!file.exists(input)) {
    fail(paste("Input file does not exist:", input))
  }

  cut_above_ground <- as.numeric(get_arg(args, "-cut", 0.25))
  min_tree_height  <- as.numeric(get_arg(args, "-height", 2))
  buffer           <- as.numeric(get_arg(args, "-buffer", 0))
  fraction         <- as.numeric(get_arg(args, "-fraction", 0.25))

  filter_str = paste("-keep_random_fraction", fraction)

  # --- Output Paths ---
  odir <- tools::file_path_sans_ext(input)
  odir <- paste0(odir, "_blender")
  if (!dir.exists(odir)) dir.create(odir)

  base <- tools::file_path_sans_ext(basename(input))
  base <- file.path(odir, base)

  out_foliage_high  <- paste0(base, "_foliage_high.laz")
  out_wood_high     <- paste0(base, "_wood_high.laz")
  out_foliage_low   <- paste0(base, "_foliage_low.laz")
  out_wood_low      <- paste0(base, "_wood_low.laz")
  out_dtm_mesh      <- paste0(base, "_dtm.obj")
  out_rings         <- paste0(base, "_rings.obj")
  its_dir           <- file.path(dirname(base), "its")

  # --- Configuration Log ---
  cat("
============ Arbor blender module =============
Input file             :", input, "
Settings
  Cut above ground (m) :", cut_above_ground, "
  Min tree height (m)  :", min_tree_height, "
  Border buffer (m)    :", buffer, "
  Filter               :", filter_str, "
====================================================
")

  # --- Processing ---
  params <- default_arbor_parameters

  cat("Reading point cloud\n")
  las <- lidR::readTLS(input, select = "xyzic", filter = filter_str)
  las <- hybrid_homogeneization(las)
  gc()

  cat("Centering on (0,0,0)\n")
  xoffset = round(mean(las$X), 1)
  yoffset = round(mean(las$Y), 1)
  zoffset = round(min(las$Z), 1)
  las@header$`X offset` = 0
  las@header$`Y offset` = 0
  las@header$`Z offset` = 0
  las$X = las$X - xoffset
  las$Y = las$Y - yoffset
  las$Z = las$Z - zoffset
  las = lidR::las_update(las)
  lidR::las_quantize(las)

  cat("Ground classification\n")
  las <- arbor_ground(las)
  gc()

  cat("DTM & height above ground\n")
  dtm <- lidR::rasterize_terrain(las, 0.1, lidR::tin())
  las <- lidR::height_above_ground(las, algorithm = dtm)
  las <- lidR::filter_poi(las, hag > cut_above_ground)
  gc()

  cat("Wood likelihood\n")
  las <- wood_likelihood(las, params)

  cat("Semantic segmentation\n")
  las <- segment_semantic(las, dtm, params)

  cat("Seeds\n")
  seeds <- find_seeds(las, params)

  cat("Instance segmentation\n")
  las <- segment_instance(las, seeds, params)

  cat("Cleaning segmentation\n")
  if (buffer > 0) las <- clip_buffer(las, -buffer)

  cat("Colorization\n")
  las <- colorize_trees(las)

  cat("Computing qsf\n")
  qsf = qsf(las)

  rings = lapply(qsf, qsm_ring)
  rings = do.call(rbind, rings)
  rings$cyl_ID = 1:nrow(rings)
  rings = set_qsm_class(rings)

  # --- Exports ---
  cat("Export qsf\n")
  qsf_write(qsf, its_dir, formats = "obj")

  cat("Export rings\n")
  qsm_write(rings, out_rings)

  cat("Export las\n")
  high = remove_small_trees(las, 2)
  low = keep_small_trees(las, 2)
  rm(las) ; gc()
  high_f = lidR::filter_poi(high, foliage > 0)
  high_w = lidR::filter_poi(high, foliage == 0)
  rm(high) ; gc()
  low_f = lidR::filter_poi(low, foliage > 0)
  low_w = lidR::filter_poi(low, foliage == 0)
  rm(low) ; gc()

  lidR::writeLAS(high_f, out_foliage_high)
  lidR::writeLAS(low_f, out_foliage_low)
  lidR::writeLAS(high_w, out_wood_high)
  lidR::writeLAS(low_w, out_wood_low)
  write_raster_to_obj(dtm, out_dtm_mesh)
}

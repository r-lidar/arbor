cmd_segment <- function(args) {

  # --- Segment Usage ---
  usage_segment <- function() {
    cat("
Usage:
  arbor segment <input.las> [options]

Mandatory:
  <input.las>             Input LAS/LAZ file

Options:
  -cut <m>                Height threshold above ground [default: 0.25]
  -height <m>             Remove trees smaller than this [default: 2]
  -buffer <m>             Buffer removed from borders [default: 5]
  -fraction <0-1>         Keep random fraction [default: 0.25]
  -tls                    Set parameters for tls
  -mls                    Set parameters for mls (default)

Export options (enabled by default):
  -no-segmented           Do not write segmented point cloud
  -no-trees               Do not write all trees
  -no-valid-trees         Do not write valid trees
  -no-dtm                 Do not write DTM raster
  -no-individual          Do not export individual trees (its/)

Other:
  -h, --help              Show this help
  -center                 Center the scene on (0,0,0)
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
  buffer           <- as.numeric(get_arg(args, "-buffer", 5))
  fraction         <- as.numeric(get_arg(args, "-fraction", 0.25))

  # Helper for export flags (local to this function)
  export_enabled <- function(name) { !paste0("-no-", name) %in% args }

  export_segmented   <- export_enabled("segmented")
  export_trees       <- export_enabled("trees")
  export_valid_trees <- export_enabled("valid-trees")
  export_dtm         <- export_enabled("dtm")
  export_individual  <- export_enabled("individual")
  center             <- has_flag(args, "-center")
  export_dtm_mesh    <- has_flag(args, "-mesh")
  tls                <- has_flag(args, "-tls")

  filter_str = paste("-keep_random_fraction", fraction)

  # --- Output Paths ---
  odir <- tools::file_path_sans_ext(input)
  odir <- paste0(odir, "_output")
  if (!dir.exists(odir)) dir.create(odir)

  base <- tools::file_path_sans_ext(basename(input))
  base <- file.path(odir, base)

  out_segmented   <- paste0(base, "_segmented.laz")
  out_trees       <- paste0(base, "_trees.laz")
  out_validtrees  <- paste0(base, "_validtrees.laz")
  out_dtm         <- paste0(base, "_dtm.tif")
  out_dtm_mesh    <- paste0(base, "_dtm.obj")
  its_dir         <- file.path(dirname(base), "its")

  # --- Configuration Log ---
  cat("
============ Arbor segmentation module =============
Input file             :", input, "
Settings
  Cut above ground (m) :", cut_above_ground, "
  Min tree height (m)  :", min_tree_height, "
  Border buffer (m)    :", buffer, "
  Filter               :", filter_str, "
Exports
  Segmented cloud      :", export_segmented, "
  All trees            :", export_trees, "
  Valid trees          :", export_valid_trees, "
  DTM                  :", export_dtm, "
  Individual trees     :", export_individual, "
====================================================
")

  # --- Processing ---
  set_lidr_threads(0)
  params <- default_arbor_parameters
  params$path_finder$max_gap = 1
  params$path_finder$k_neighborhood_connectivity = 20

  if (tls)

  cat("Reading point cloud\n")
  las <- readTLS(input, select = "xyzic", filter = filter_str)
  las <- hybrid_homogeneization(las)
  gc()

  if (center)
  {
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
    las = las_update(las)
    las_quantize(las)
  }

  cat("Ground classification\n")
  las <- arbor_ground(las)
  gc()

  cat("DTM & height above ground\n")
  dtm <- rasterize_terrain(las, 0.1, tin())
  las <- height_above_ground(las, algorithm = dtm)
  las <- filter_poi(las, hag > cut_above_ground)
  gc()

  cat("Wood likelihood\n")
  las <- wood_likelihood(las, params)

  cat("Semantic segmentation\n")
  las  <- segment_semantic(las, dtm, params)

  cat("Seeds\n")
  seeds <- find_seeds(las, params)

  cat("Instance segmentation\n")
  las <- segment_instance(las, seeds, params)

  cat("Cleaning segmentation\n")
  trees <- remove_small_trees(las, max_height = min_tree_height)
  trees <- fix_small_isolated_low_clusters(trees)
  trees <- clip_buffer(trees, -buffer)

  cat("Colorization\n")
  las         <- colorize_trees(las)
  trees       <- colorize_trees(trees)

  # --- Exports ---
  cat("Exports\n")
  if (export_segmented)   writeLAS(las, out_segmented)
  if (export_trees)       writeLAS(trees, out_trees)
  if (export_dtm)         writeRaster(dtm, out_dtm, overwrite = TRUE)
  if (export_dtm_mesh)    arbor:::write_raster_to_obj(dtm, out_dtm_mesh)

  name = tools::file_path_sans_ext(basename(input))

  if (export_individual) {
    dir.create(its_dir, showWarnings = FALSE, recursive = TRUE)
    for (i in unique(trees$treeID)) {
      tree <- filter_poi(trees, treeID == i)
      writeLAS(tree, file.path(its_dir, paste0(name, "_tree_", i, ".las")))
    }
  }
}

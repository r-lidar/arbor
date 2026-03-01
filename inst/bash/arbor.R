#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lidR)
  library(arbor)
  library(terra)
  library(tools)
})

# ------------------------------------------------------------
# Global Helpers
# ------------------------------------------------------------

# Parse a value associated with a flag (e.g., --buffer 5)
get_arg <- function(args, flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 1 && idx < length(args)) {
    return(args[idx + 1])
  }
  default
}

# Check if a flag exists (boolean)
has_flag <- function(args, flag) {
  flag %in% args
}

# Fail safely with a message
fail <- function(msg) {
  cat("Error:", msg, "\n", file = stderr())
  quit(status = 1)
}

# Main Usage Help
usage_main <- function() {
  cat("
Usage:
  arbor <command> [arguments]

Commands:
  segment      Segment a point cloud (LAS/LAZ)
  qsm          Run QSM on a folder or file
  report       Produce a pdf report

Options:
  -h, --help   Show help for a specific command

Examples:
  arbor segment plot.laz --no-dtm
  arbor qsm ./my_trees/ -ncores 4
")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------
# Command: Segment
# ------------------------------------------------------------

run_segment <- function(args) {

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

Export options (enabled by default):
  -no-segmented           Do not write segmented point cloud
  -no-trees               Do not write all trees
  -no-valid-trees         Do not write valid trees
  -no-dtm                 Do not write DTM raster
  -no-individual          Do not export individual trees (its/)

Other:
  -h, --help               Show this help
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
  classifier = csf(rigidness = 1, class_threshold = 0.05, cloth_resolution = 0.1)
  las <- classify_ground(las, classifier, last_returns = FALSE)
  ground <- filter_poi(las, Classification == LASGROUND)
  ground <- decimate_points(ground, lowest(0.25))
  ground <- classify_noise(ground, sor(k = 10, m = 2))
  ground <- remove_noise(ground)
  ground$Classification <- LASGROUND
  gc()

  cat("DTM & height above ground\n")
  dtm <- rasterize_terrain(ground, 0.1, tin())
  las <- height_above_ground(las, algorithm = dtm)
  las <- filter_poi(las, hag > cut_above_ground)
  gc()

  cat("Wood likelihood\n")
  las    <- wood_likelihood(las, params)

  cat("Semantic segmentation\n")
  las    <- segment_semantic(las, dtm, params)

  cat("Seeds\n")
  seeds  <- find_seeds(las, params)

  cat("Instance segmentation\n")
  las    <- segment_instance(las, seeds, params)

  cat("Cleaning segmentation\n")
  trees <- remove_small_trees(las, max_height = min_tree_height)
  trees <- fix_small_isolated_low_clusters(trees)
  valid_trees <- clip_buffer(trees, -buffer)


  cat("Colorization\n")
  las         <- colorize_trees(las)
  trees       <- colorize_trees(trees)
  valid_trees <- colorize_trees(valid_trees)

  # --- Exports ---
  cat("Exports\n")
  if (export_segmented)   writeLAS(las, out_segmented)
  if (export_trees)       writeLAS(trees, out_trees)
  if (export_valid_trees) writeLAS(valid_trees, out_validtrees)
  if (export_dtm)         writeRaster(dtm, out_dtm, overwrite = TRUE)
  if (export_dtm_mesh)    arbor:::write_raster_to_obj(dtm, out_dtm_mesh)

  name = tools::file_path_sans_ext(basename(input))

  if (export_individual) {
    dir.create(its_dir, showWarnings = FALSE, recursive = TRUE)
    for (i in unique(valid_trees$treeID)) {
      tree <- filter_poi(valid_trees, treeID == i)
      writeLAS(tree, file.path(its_dir, paste0(name, "_tree_", i, ".las")))
    }
  }
}

# ------------------------------------------------------------
# Command: QSM
# ------------------------------------------------------------

run_qsm <- function(args) {

  # --- QSM Usage ---
  usage_qsm <- function() {
    cat("
Usage:
  arbor qsm <input> [options]

Mandatory:
  <input>           Input folder or file

Options:
  -o, --output <dir> Output directory [default: input directory]
  -ncores <int>      Number of cores [default: half of available]
  -overwrite         Overwrite existing files
  -csv               Export CSV [default: on if nothing provided]
  -obj               Export OBJ [default: on if nothing provided]
  -ply               Export PLY

Example:
  arbor qsm ./its/ -csv -obj --ncores 4
")
    quit(save = "no", status = 0)
  }

  if (length(args) == 0 || has_flag(args, "-h") || has_flag(args, "--help")) {
    usage_qsm()
  }

  # --- Parse Inputs ---
  # Prioritize positional argument for input
  ifiles <- NULL
  if (!startsWith(args[1], "-")) {
    ifiles <- args[1]
  } else {
    # Fallback to flags if user uses legacy style
    ifiles <- get_arg(args, "-i")
    if (is.null(ifiles)) ifiles <- get_arg(args, "-input")
  }

  if (is.null(ifiles)) fail("Missing input file or folder")
  if (!file.exists(ifiles)) fail(paste("Input does not exist:", ifiles))
  ifiles <- normalizePath(ifiles, mustWork = TRUE)

  # Output dir
  odir <- file.path(dirname(ifiles), "qsm")

  # Override output if flag present
  out_flag <- get_arg(args, "-o")
  if (is.null(out_flag)) out_flag <- get_arg(args, "--output")
  if (!is.null(out_flag)) odir <- normalizePath(out_flag, mustWork = FALSE)

  # Formats
  formats <- character()
  if (has_flag(args, "-csv")) formats <- c(formats, "csv")
  if (has_flag(args, "-obj")) formats <- c(formats, "obj")
  if (has_flag(args, "-ply")) formats <- c(formats, "ply")
  if (length(formats) == 0) formats <- c("csv", "obj") # Default

  overwrite <- has_flag(args, "-overwrite")

  # Cores
  ncores_arg <- get_arg(args, "-ncores")
  if (is.null(ncores_arg)) ncores_arg <- get_arg(args, "-ncores")

  if (!is.null(ncores_arg)) {
    ncores <- as.integer(ncores_arg)
  } else {
    ncores <- floor(parallel::detectCores() / 2)
  }

  if (is.na(ncores) || ncores < 1) fail("--ncores must be a positive integer")

  cat("
================= Arbor QSM module =================
Input       :", paste(ifiles, collapse = ", "), "
Output      :", odir, "
Settings
  Formats   :", paste(formats, collapse = ", "), "
  Cores     :", ncores, "
=====================================================
")

  res <- qsf(
    input     = ifiles,
    odir      = odir,
    formats   = formats,
    overwrite = overwrite,
    ncores    = ncores
  )

  log <- qsf_log(res)

  #for (i in which(has_msg))
  #{
  #  if (isTRUE(res$success[i])) {
  #    cat("\033[33m", "WARNING | ", res$name[i], " | ", res$message[i], "\033[0m\n", sep = "")
  #  } else {
  #    cat("\033[31m", "ERROR   | ", res$name[i], " | ", res$message[i], "\033[0m\n", sep = "")
  #  }
  #}
}

# ----------
# Report
# ----------

run_report = function(args)
{
  usage_report = function()
  {
    cat("
Usage:
  arbor report <input_dir> <output_pdf>
")
    quit(save = "no", status = 0)
  }

  if (length(args) < 2 || has_flag(args, "-h") || has_flag(args, "--help")) {
    usage_report()
  }

  input_dir <- normalizePath(args[1])
  output_pdf <- normalizePath(args[2], mustWork = FALSE)

  if (!dir.exists(input_dir)) {
    stop("Input directory does not exist: ", input_dir)
  }

  rmarkdown::render(
    input  = system.file("bash", "report_template.Rmd", package="arbor"),
    output_file = basename(output_pdf), # only filename
    output_dir = dirname(output_pdf), # directory to write to
    params = list(
      idir = input_dir
    ),
    quiet = FALSE
  )
}

# ------------------------------------------------------------
# Main Dispatcher
# ------------------------------------------------------------

all_args <- commandArgs(trailingOnly = TRUE)

if (length(all_args) == 0) {
  usage_main()
}

command  <- all_args[1]
sub_args <- all_args[-1]

t0 = arbor:::tic()

if (command == "segment") {
  run_segment(sub_args)
} else if (command == "qsm") {
  run_qsm(sub_args)
} else if (command == "report") {
  run_report(sub_args)
} else if (command %in% c("-h", "--help")) {
  usage_main()
} else {
  cat(paste0("Unknown command: '", command, "'\n"))
  usage_main()
}

arbor:::toc(t0, space = "")

# @file cmd_segment.R
# Project: Arbor
#
# Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Options known by `arbor segment`
SEGMENT_VALUE_FLAGS <- c("-cut", "-height", "-buffer", "-fraction", "-set_epsg", "-aoi")
SEGMENT_BOOL_FLAGS  <- c("-no-segmented", "-no-trees", "-no-dtm", "-center", "-mesh", "-memory", "-h", "--help")

# Read an area of interest: a single polygon without hole
read_aoi <- function(file)
{
  if (!file.exists(file)) stop("AOI file does not exist: ", file, call. = FALSE)
  aoi <- sf::st_zm(sf::st_geometry(sf::st_read(file, quiet = TRUE)))

  if (length(aoi) == 1L && inherits(aoi[[1]], "MULTIPOLYGON") && length(aoi[[1]]) == 1L)
    aoi <- sf::st_cast(aoi, "POLYGON")

  if (length(aoi) != 1L || !inherits(aoi[[1]], "POLYGON") || length(aoi[[1]]) != 1L)
    stop("The AOI must be a single polygon without hole (multipolygons are not supported).", call. = FALSE)

  aoi
}

cmd_segment <- function(args, strict = TRUE) {

  treeID <- UserData <- NULL

  # --- Segment Usage ---
  usage_segment <- function() {
    cat("
Usage:
  arbor segment <input.las> [options]

Mandatory:
  <input.las>             Input LAS/LAZ file

Options:
  -cut <m>                Points below this height above ground are ignored [default: 0.25]
  -height <m>             Remove trees smaller than this [default: 2]
  -buffer <m>             Width of the border strip whose trees are flagged [default: 10]
  -fraction <0-1>         Fraction of points randomly kept at reading [default: 0.25]
  -set_epsg 1234          EPSG code to georeference the dataset if missing
  -aoi file.shp           GDAL readable polygon (single polygon, no hole) defining the area
                          of interest. If missing, the convex hull of the point cloud is used

Export options (enabled by default):
  -no-segmented           Do not write segmented point cloud
  -no-trees               Do not write all trees
  -no-dtm                 Do not write DTM raster

Other:
  -h, --help              Show this help
  -center                 Center the scene on (0,0,0) (the AOI is shifted accordingly)
  -mesh                   Export DTM as mesh (.obj)
  -memory                 Monitor memory usage and write it to a text file
")
    quit(save = "no", status = 0)
  }

  if (length(args) == 0 || has_flag(args, "-h") || has_flag(args, "--help")) {
    usage_segment()
  }

  if (strict) check_args(args, SEGMENT_VALUE_FLAGS, SEGMENT_BOOL_FLAGS)

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
  buffer           <- as.numeric(get_arg(args, "-buffer", 10))
  fraction         <- as.numeric(get_arg(args, "-fraction", 0.25))
  epsg             <- as.numeric(get_arg(args, "-set_epsg", 0))
  faoi             <- get_arg(args, "-aoi", NA_character_)

  if (is.na(buffer) || buffer < 0) fail("-buffer must be a positive number.")
  if (is.na(fraction) || fraction <= 0 || fraction > 1) fail("-fraction must be in ]0, 1].")
  if (is.na(epsg)) fail("-set_epsg must be an integer EPSG code.")

  # numeric (negative buffer from the convex hull) or sf polygon (area of interest)
  aoi <- -buffer
  if (!is.na(faoi))
  {
    aoi <- tryCatch(read_aoi(faoi), error = function(e) fail(conditionMessage(e)))
  }

  # Helper for export flags (local to this function)
  export_enabled <- function(name) { !paste0("-no-", name) %in% args }

  export_segmented   <- export_enabled("segmented")
  export_trees       <- export_enabled("trees")
  export_dtm         <- export_enabled("dtm")

  center             <- has_flag(args, "-center")
  export_dtm_mesh    <- has_flag(args, "-mesh")
  track_memory       <- has_flag(args, "-memory")


  filter_str = paste("-keep_random_fraction", fraction)

  # --- Output Paths ---
  odir <- tools::file_path_sans_ext(input)
  odir <- paste0(odir, "_output")
  if (!dir.exists(odir)) dir.create(odir)

  base <- tools::file_path_sans_ext(basename(input))
  base <- file.path(odir, base)

  out_segmented   <- paste0(base, "_segmented.laz")
  out_trees       <- paste0(base, "_trees.laz")
  out_dtm         <- paste0(base, "_dtm.tif")
  out_dtm_mesh    <- paste0(base, "_dtm.obj")
  out_memory      <- paste0(base, "_memory.txt")
  its_dir         <- file.path(dirname(base), "its")

  crs = sf::NA_crs_
  if (epsg != 0) crs <- sf::st_crs(paste0("EPSG:", epsg))

  if (is.numeric(aoi))
    aoistr = "convex hull"
  else
    aoistr = faoi

  # --- Configuration Log ---
  cat("
============ Arbor segmentation module =============
Input file             :", input, "
Settings
  Cut above ground (m) :", cut_above_ground, "
  Min tree height (m)  :", min_tree_height, "
  Buffer (m)           :", buffer, "
  Area of interest     :", aoistr, "
  Filter               :", filter_str, "
  Centering            :", center, "
Exports
  Segmented cloud      :", export_segmented, "
  All trees            :", export_trees, "
  DTM                  :", export_dtm, "
====================================================
")

  if (track_memory) f <- mem_monitor(out_memory, 1)

  # --- Processing ---
  params <- arbor_parameters_default
  params$global$cut_above_ground <- cut_above_ground

  if (track_memory) phase("Read", f)
  cat("Reading point cloud\n")
  las <- read_tls(input, select = "xyzic", filter = filter_str)
  gc()

  if (!is.na(crs))
  {
    sf::st_crs(las) <- crs
  }

  if (!is.numeric(aoi))
  {
    if (is.na(sf::st_crs(aoi)) || is.na(sf::st_crs(las)))
      sf::st_crs(aoi) <- sf::st_crs(las)
    else if (sf::st_crs(aoi) != sf::st_crs(las))
      aoi <- sf::st_transform(aoi, sf::st_crs(las))

    cat("Clipping the point cloud\n")
    bpol <- sf::st_buffer(aoi, buffer)
    if (track_memory) phase("Clip", f)
    las <- lidR::clip_roi(las, bpol)
  }

  if (track_memory) phase("Sampling", f)
  cat("Hybrid homogeneization\n")
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
    las = lidR::las_update(las)
    lidR::las_quantize(las)

    # The AOI must follow the point cloud
    if (!is.numeric(aoi))
    {
      crs_aoi <- sf::st_crs(aoi)
      aoi <- sf::st_set_crs(sf::st_geometry(aoi) - c(xoffset, yoffset), crs_aoi)
    }
  }

  if (track_memory) phase("Ground", f)
  cat("Ground classification & height above ground\n")
  las <- segment_ground(las)
  gc()

  if (track_memory) phase("DTM", f)
  cat("DTM\n")
  dtm <- lidR::rasterize_terrain(las, 0.1, lidR::tin())
  gc()

  if (track_memory) phase("Likelihood", f)
  cat("Wood likelihood\n")
  las <- wood_likelihood(las, params)

  if (track_memory) phase("Semantic", f)
  cat("Semantic segmentation\n")
  las  <- segment_semantic(las, params)

  if (track_memory) phase("Seeds", f)
  cat("Seeds\n")
  see <- find_seeds(las, params)

  if (track_memory) phase("Instance", f)
  cat("Instance segmentation\n")
  las <- segment_instance(las, see, params)

  if (track_memory) phase("Flag", f)
  cat("Cleaning segmentation\n")
  las <- flag_small_trees(las, max_height = min_tree_height)
  las <- flag_buffer(las, see, aoi)

  if (track_memory) phase("RGB", f)
  cat("Colorization\n")
  las <- colorize_trees(las)

  trees <- lidR::filter_poi(las, UserData == ARBORTREE)

  if (track_memory) mem_monitor_stop(f)

  # --- Exports ---
  cat("Exports\n")
  if (export_segmented)   lidR::writeLAS(las, out_segmented)
  if (export_trees)       lidR::writeLAS(trees, out_trees)
  if (export_dtm)         terra::writeRaster(dtm, out_dtm, overwrite = TRUE)
  if (export_dtm_mesh)    write_raster_to_obj(dtm, out_dtm_mesh)

  invisible(list(
    dir        = odir,
    segmented  = if (export_segmented) out_segmented else NA_character_,
    trees      = if (export_trees)     out_trees     else NA_character_,
    dtm        = if (export_dtm)       out_dtm       else NA_character_,
    dtm_mesh   = if (export_dtm_mesh)  out_dtm_mesh  else NA_character_
  ))
}

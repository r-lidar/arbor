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

cmd_segment <- function(args) {

  treeID <- UserData <- NULL

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
====================================================
")

  # --- Processing ---
  params <- arbor_parameters_default
  params$path_finder$max_gap = 1
  params$path_finder$k_neighborhood_connectivity = 20

  if (tls)
  {}

  cat("Reading point cloud\n")
  las <- lidR::readTLS(input, select = "xyzic", filter = filter_str)
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
  }

  cat("Ground classification & height above ground\n")
  las <- segment_ground(las)
  gc()

  cat("DTM\n")
  dtm <- lidR::rasterize_terrain(las, 0.1, lidR::tin())
  gc()

  cat("Wood likelihood\n")
  las <- wood_likelihood(las, params)

  cat("Semantic segmentation\n")
  las  <- segment_semantic(las, params)

  cat("Seeds\n")
  see <- find_seeds(las, params)

  cat("Instance segmentation\n")
  las <- segment_instance(las, see, params)

  cat("Cleaning segmentation\n")
  las <- flag_small_trees(las, max_height = min_tree_height)
  las <- flag_buffer(trees, las, -buffer)

  cat("Colorization\n")
  las <- colorize_trees(las)

  trees <- lidR::filter_poi(las, UserData == ARBORTREE)

  # --- Exports ---
  cat("Exports\n")
  if (export_segmented)   lidR::writeLAS(las, out_segmented)
  if (export_trees)       lidR::writeLAS(trees, out_trees)
  if (export_dtm)         terra::writeRaster(dtm, out_dtm, overwrite = TRUE)
  if (export_dtm_mesh)    write_raster_to_obj(dtm, out_dtm_mesh)
}

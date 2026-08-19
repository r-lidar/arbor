# @file cmd_blender.R
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

cmd_blender <- function(args) {

  UserData <- NULL

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
  -h, --help              Show this help
  -center                 Center the scene on (0,0,0)
")
    quit(save = "no", status = 0)
  }

  if (length(args) == 0 || has_flag(args, "-h") || has_flag(args, "--help")) {
    usage_blender()
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
  params <- arbor_parameters_default

  cat("Reading point cloud\n")
  las <- lidR::readTLS(input, select = "xyzic", filter = filter_str)
  t0  <- Sys.time()
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

  cat("Ground classification & height above ground\n")
  las <- segment_ground(las)
  gc()

  cat("DTM\n")
  dtm <- lidR::rasterize_terrain(las, 0.1, lidR::tin())
  gc()

  cat("Wood likelihood\n")
  las <- wood_likelihood(las, params)
  gc()

  cat("Semantic segmentation\n")
  las <- segment_semantic(las, params)
  gc()

  cat("Seeds\n")
  seeds <- find_seeds(las, params)
  gc()

  cat("Instance segmentation\n")
  las <- segment_instance(las, seeds, params)
  gc()

  cat("Cleaning segmentation\n")
  if (buffer > 0)
  {
    treeID <- NULL
    las <- flag_buffer(las, seeds, -buffer)
    las <- lidR::filter_poi(las, UserData != ARBORBUFFER)
  }

  cat("Colorization\n")
  las <- colorize_trees(las)
  gc()

  cat("Computing qsf\n")
  qsf = qsf(las)

  tf = Sys.time()
  dt = difftime(tf, t0)
  record_entry(las, dt, input)

  rings = lapply(qsf, qsm_ring)
  rings = do.call(rbind, rings)
  rings$cyl_ID = 1:nrow(rings)
  rings = as_qsm(rings)

  # --- Exports ---
  cat("Export qsf\n")
  qsf_write(qsf, its_dir, formats = "obj")

  cat("Export rings\n")
  write_cylinders_to_obj(rings, out_rings)

  cat("Export las\n")
  foliage <- NULL
  high = flag_small_trees(las, 2)
  high = lidR::filter_poi(las, UserData == ARBORTREE)
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

write_cylinders_to_obj <- function(df, file_dest, sides = 12) {
  # Internal helper for cross products
  cross_prod <- function(a, b) {
    c(a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1])
  }

  all_output <- character(nrow(df) * 2 + 2)
  all_output[1] <- "# Exported Tree Cylinders"

  v_list <- list()
  f_list <- list()
  v_offset <- 1

  for (i in seq_len(nrow(df))) {
    # Extract points and radius
    p1 <- as.numeric(df[i, c("startX", "startY", "startZ")])
    p2 <- as.numeric(df[i, c("endX", "endY", "endZ")])
    r  <- df$radius[i]

    # 1. Calculate orientation
    axis <- p2 - p1
    len <- sqrt(sum(axis^2))
    if(len == 0) next # Skip degenerate cylinders
    axis <- axis / len

    # 2. Create orthogonal basis for the cylinder circles
    # Find a vector not parallel to the axis
    tmp <- if (abs(axis[1]) < 0.9) c(1, 0, 0) else c(0, 1, 0)
    u <- cross_prod(axis, tmp)
    u <- u / sqrt(sum(u^2))
    v <- cross_prod(axis, u)

    # 3. Generate vertices for both ends
    verts <- matrix(0, nrow = sides * 2, ncol = 3)
    for (s in 0:(sides - 1)) {
      ang <- 2 * pi * s / sides
      cp <- (cos(ang) * u + sin(ang) * v) * r
      verts[s + 1, ] <- p1 + cp
      verts[s + 1 + sides, ] <- p2 + cp
    }

    # Format vertices as OBJ strings
    v_list[[i]] <- apply(verts, 1, function(row) sprintf("v %.6f %.6f %.6f", row[1], row[2], row[3]))

    # 4. Generate faces (connecting the two rings)
    faces <- character(sides)
    for (s in 1:sides) {
      v1 <- s + v_offset - 1
      v2 <- (s %% sides + 1) + v_offset - 1
      v3 <- v1 + sides
      v4 <- v2 + sides
      # Define two triangles for each quad segment
      faces[s] <- sprintf("f %d %d %d\nf %d %d %d", v1, v2, v3, v2, v4, v3)
    }
    f_list[[i]] <- faces

    # Update global vertex offset
    v_offset <- v_offset + (sides * 2)
  }

  # Combine and write to disk
  final_content <- c("# Vertices", unlist(v_list), "# Faces", unlist(f_list))
  writeLines(final_content, file_dest)
  message(paste("Successfully wrote OBJ to:", file_dest))
}

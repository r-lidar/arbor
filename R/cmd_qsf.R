# @file cmd_qsm.R
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

# Options known by `arbor qsf`
QSF_VALUE_FLAGS <- c("-o", "--output", "-allometry", "-set_epsg", "-i", "-input")
QSF_BOOL_FLAGS  <- c("-csv", "-obj", "-ply", "-h", "--help")

cmd_qsf <- function(args, strict = TRUE) {

  # --- QSM Usage ---
  usage_qsm <- function() {
    cat("
Usage:
  arbor qsf <input> [options]

Mandatory:
  <input>           Input LAS file segmented with arbor

Options:
  -o, --output <dir>  Output directory [default: input directory]
  -allometry <name>   Choose an allometry name
  -set_epsg 1234      EPSG code to georeference the dataset if missing
  -csv                Export CSV
  -obj                Export OBJ
  -ply                Export PLY

The <name>.qsf manifest and the qsm/ folder are always written (existing files are overwritten).

Example:
  arbor qsf segmented.laz
")
    quit(save = "no", status = 0)
  }

  if (length(args) == 0 || has_flag(args, "-h") || has_flag(args, "--help")) {
    usage_qsm()
  }

  if (strict) check_args(args, QSF_VALUE_FLAGS, QSF_BOOL_FLAGS)

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
  odir <- dirname(ifiles)

  # Override output if flag present
  out_flag <- get_arg(args, "-o")
  if (is.null(out_flag)) out_flag <- get_arg(args, "--output")
  if (!is.null(out_flag)) odir <- normalizePath(out_flag, mustWork = FALSE)

  epsg <- as.numeric(get_arg(args, "-set_epsg", 0))
  crs = sf::NA_crs_
  if (epsg != 0) crs <- sf::st_crs(paste0("EPSG:", epsg))

  # Parameters
  allometry = get_arg(args, "-allometry", NULL)
  params = arbor_parameters_default
  if (!is.null(allometry)) params$qsm$allometry_name = allometry

  # Formats. The .qsf manifest and the .qsm files are always written (once), the other
  # formats are optional.
  formats <- character()
  if (has_flag(args, "-csv")) formats <- c(formats, "csv")
  if (has_flag(args, "-obj")) formats <- c(formats, "obj")
  if (has_flag(args, "-ply")) formats <- c(formats, "ply")

  cat("
================= Arbor QSM module =================
Input       :", paste(ifiles, collapse = ", "), "
Output      :", odir, "
Settings
  Formats   :", paste(c("qsm", formats), collapse = ", "), "
  Allometry :", params$qsm$allometry_name, "
=====================================================
")

  cat("Reading point cloud\n")
  las <- lidR::readLAS(ifiles)

  cat("Computing QSMs\n")
  res <- qsf(las, params = params)

  if (!is.na(crs))
    sf::st_crs(res) <- crs

  cat("Writing QSMs\n")
  if (!dir.exists(odir)) dir.create(odir, recursive = TRUE)
  qsf_file <- file.path(odir, paste0(tools::file_path_sans_ext(basename(ifiles)), ".qsf"))

  # Single call writing odir/<name>.qsf and odir/qsm/*.qsm
  qsf_write(res, qsf_file)

  # Optional extra formats: odir/<format>/*.<format>
  if (length(formats) > 0)
    qsf_write(res, odir, formats)

  invisible(list(
    dir      = odir,
    qsf_file = qsf_file
  ))
}

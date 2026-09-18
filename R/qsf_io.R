# @file qsf_io.R
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

#' Write a QSF to Files
#'
#' Writes a Quantitative Structure Forest (QSF) to files in either OBJ, PLY or CSV format,
#' based on the file extension. By default, each QSM is written in its own file inside `dir`.
#' Alternatively, if `dir` is given as a single file path whose extension matches `formats`
#' (e.g. `"forest.obj"` or `"forest.ply"`), all the QSM are merged and written into that
#' single OBJ or PLY file. This single-file mode is not supported for the `.qsm` format.
#' In single-file mode, each QSM is kept as an individually named object: OBJ files get one
#' `o <name>` group per QSM, while PLY files get an extra per-face `object_id` property (and a
#' `comment object <id> <name>` header line) identifying which QSM each face belongs to.
#' OBJ has no binary variant and is always written as ASCII regardless of `binary`; PLY supports
#' both ASCII (`binary = FALSE`) and binary little-endian (`binary = TRUE`).
#'
#' Supported formats:
#' * `.qsm` native binary format
#' * `.ply` or `.obj` or `.stl`: writes the QSM as a mesh file
#' * `.csv` or `.txt`: writes the QSM as a ASCII table
#'
#' @param qsf A QSF object to be written.
#' @param dir A string giving the directory to the output files, or a single file path
#' (with an extension matching `formats`) to write all the QSM into one OBJ or PLY file.
#' @param formats the format (e.g. "qsm", "ply", "obj", "csv", "txt", "stl").
#' @param binary Boolean. Used if the format supports ASCII or binary
#'
#' @export
#' @export
#' @md
qsf_write = function(qsf, dir, formats = c("qsm", "obj"), binary = TRUE)
{
  dir = normalizePath(dir, mustWork = FALSE)
  for (format in formats) qsf_write_cpp(qsf, dir, format, binary)
  return(invisible(TRUE))
}

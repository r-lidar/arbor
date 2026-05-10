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
#' based on the file extension. Each QSM is written in its own file.
#'
#' Supported formats:
#' * `.ply` or `.obj` or `.stl`: writes the QSM as a mesh file
#' * `.csv` or `.txt`: writes the QSM as a table
#'
#' @param qsf A QSF object to be written.
#' @param dir A string giving the director to the output files.
#' @param formats the format (e.g., ".ply", "obj", ".csv", ".txt", ".stl").
#' @param binary Boolean. Used if the format supports ASCII or binary
#'
#' @export
#' @export
#' @md
qsf_write = function(qsf, dir, formats = c("csv", "obj"), binary = TRUE)
{
  dir = normalizePath(dir, mustWork = FALSE)
  for (format in formats) qsf_write_cpp(qsf, dir, format, binary)
  return(invisible(TRUE))
}

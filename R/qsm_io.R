# @file qsm_io.R
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

#' Write a QSM to File
#'
#' Writes a Quantitative Structure Model (QSM) to a file in either OBJ, PLY or CSV format,
#' based on the file extension.
#'
#' Supported formats:
#' * `.ply` or `.obj` or `.stl`: writes the QSM as a mesh file
#' * `.csv` or `.txt`: writes the QSM as a table
#'
#' @param qsm A QSM object to be written.
#' @param file A string giving the path to the output file. The file extension determines
#'   the format (e.g., ".ply", "obj", ".csv", ".txt", ".stl").
#' @param binary Boolean. Used if the format supports ASCII or binary
#'
#' @export
#'
#' @examples
#' \dontrun{
#' qsm_write(qsm, "tree.ply")
#' qsm_write(qsm, "tree.csv")
#' }
#' @export
#' @md
qsm_write = function(qsm, file, binary = TRUE)
{
  file = normalizePath(file, mustWork = FALSE)
  qsm_write_cpp(qsm, file, binary)
  return(invisible(TRUE))
}

#' Read QSM Data from File
#'
#' Loads QSM segment data from a CSV file and unify naming convention. So the function reads Computree
#' TreeQSM, arbor in the same format.
#'
#' @param x A character string specifying the file path.
#' @return A data frame of QSM segment data.
#' @export
#' @md
qsm_read = function(x)
{
  qsm = data.table::fread(x)
  qsm = unify_names(qsm)
  name = tools::file_path_sans_ext(basename(x))
  attr(qsm, "ID") = name
  qsm = as_qsm(qsm)
  return(qsm)
}

unify_names <- function(qsm)
{
  # Define the renaming mapping: old = new
  name_map <- c(
    "branchOrder" = "branch_order",
    "parent_segment_ID" = "parent_ID",
    "segment_ID" = "cyl_ID"
  )

  # Replace column names only if they exist in the data.frame
  current_names <- names(qsm)
  new_names <- ifelse(current_names %in% names(name_map), name_map[current_names], current_names)
  names(qsm) <- new_names
  qsm
}


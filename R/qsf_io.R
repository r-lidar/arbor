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
#' Writes a Quantitative Structure Forest (QSF), i.e. a collection of QSMs
#' to files in QSM, OBJ, PLY, STL, CSV, or TXT format.
#'
#' The function operates in two distinct modes:
#' \itemize{
#'   \item \strong{Directory mode (recommended):} If `path` is a directory (or has
#'     no file extension), each QSM in the forest is written into its own separate
#'     file inside `path/format/`, named after the QSM's name (or its id if it has
#'     no name). Multiple formats can be exported at once by passing several values
#'     in `format`. Every tree is always fully separate and identifiable in this
#'     mode, for every format, because each tree gets its own file.
#'   \item \strong{Single-file mode:} If `path` points to a specific file path
#'     (i.e., its extension is non-empty), a single file is written containing all
#'     QSMs together. The file format is inferred strictly from the file
#'     extension of `path`, and the `format` argument is ignored. This is only
#'     supported for the \code{obj}, \code{ply}, and \code{stl} mesh formats.
#'     \code{qsm}, \code{csv}, and \code{txt} do not support combining several
#'     trees into one file and will raise an error if used this way. Even among the
#'     supported formats, how well individual trees stay distinguishable inside the
#'     combined file varies a lot (see "Limits of single-file mode" below), so this
#'     mode is best reserved for cases where you specifically need one file to hand
#'     to some other piece of software (e.g. importing a whole forest into Blender
#'     in one go).
#' }
#'
#' Supported formats:
#' \itemize{
#'   \item \code{.qsm}: Native binary format. Directory mode only. Also produces
#'     a sibling \code{.qsf} manifest file that indexes every \code{.qsm} file
#'     written, by relative path (see \link{qsf_read}).
#'   \item \code{.obj}, \code{.ply}, or \code{.stl}: Mesh formats. Available in
#'     both directory and single-file mode.
#'   \item \code{.csv} or \code{.txt}: ASCII table formats. Directory mode only.
#' }
#'
#' @section Limits of single-file mode:
#'
#' When several QSMs are combined into one \code{obj}, \code{ply}, or \code{stl}
#' file, each format has a different (and sometimes very limited) way of keeping
#' individual trees identifiable:
#' \itemize{
#'   \item \strong{OBJ:} Full support. Each QSM becomes its own named `o` object
#'     inside the file. This is NOT respected by all OBJ readers. E.g. Blender
#'     respects the OBJ specifications and load individual trees but not
#'     CloudCompare that loads one massive undifferentiated mesh.
#'   \item \strong{PLY:} Partial support. All QSMs share one vertex/face list (a
#'     single mesh), but every face carries a `qsm_id` integer property, and the
#'     header records one comment per QSM with its id and name. Software that lets
#'     you inspect or filter by custom face properties can recover the original
#'     trees; a generic PLY viewer will just show one undifferentiated mesh.
#'   \item \strong{STL (ASCII, \code{binary = FALSE}):} Weak, non-standard support.
#'     Each QSM is written as its own `solid <name> ... endsolid <name>` block.
#'     This is a common convention, not part of the STL standard, so support is
#'     inconsistent: some tools import every block as a separate named object,
#'     others only read the first block, and others merge everything into one
#'     mesh. Always test with your target software before relying on this.
#'   \item \strong{STL (binary, \code{binary = TRUE}, the default):} No support at
#'     all. Binary STL has no concept of an object name or grouping - it is just a
#'     flat, unlabelled list of triangles. Combining multiple QSMs into one binary
#'     STL file merges them into a single anonymous mesh with no way to tell which
#'     triangle came from which tree. Use \code{binary = FALSE} (ASCII STL), or
#'     switch to \code{obj} or \code{ply}, if keeping trees separate matters.
#' }
#'
#' @param qsf A QSF object to be written.
#' @param path A string giving the output directory path (directory mode), or a
#'   file path ending in \code{.obj}, \code{.ply}, or \code{.stl} for single-file
#'   mode (see Details).
#' @param format A character vector specifying the export format(s) (e.g., \code{"qsm"},
#'   \code{"obj"}, \code{"ply"}, \code{"stl"}, \code{"csv"}, \code{"txt"}). Defaults to
#'   \code{c("qsm", "obj")}. Ignored when \code{path} specifies a single-file output.
#' @param binary Logical. Indicates whether to output binary or ASCII where supported.
#'   Defaults to \code{TRUE}. For \code{stl} in particular this also controls whether
#'   trees stay distinguishable in single-file mode - see Details.
#'
#' @return Returns \code{TRUE} invisibly upon completion.
#' @seealso \link{qsf_read}
#' @examples
#' \dontrun{
#' # Single-file mode: maybe the simplest use case. Writes forest.qsf file that references 
#' # the .qsm files in the subdirectory "qsm/". In this case it does produces numerous
#' # files because .qsf file format only contains metadata.
#' qsf_write(qsf, "forest.qsf")
#' 
#' # Directory mode: one file per tree.
#' # Writes forest/qsm/<name>.qsm and forest/obj/<name>.obj for every tree.
#' # Writes also path/to/directory/directory.qsf to reference all the .qsm file
#' qsf_write(qsf, "path/to/directory/", format = c("qsm", "obj"))
#'
#' # Directory mode: several formats at once, ASCII where applicable.
#' qsf_write(qsf, "forest", format = c("ply", "csv"), binary = FALSE)
#'
#' # Single-file mode: with OBJ format is the safest choice. Write a single OBJ file. Every 
#' # tree keeps its name as a separate 'o' object. Be careful that  NOT all software will 
#' # import as such. Blender supports it not CloudCompare.
#' qsf_write(qsf, "forest_all_trees.obj")
#'
#' # Single-file mode: with PLY format. One mesh, but trees can still be told apart via the
#' # per-face 'qsm_id' property if your software can filter on it.
#' qsf_write(qsf, "forest_all_trees.ply")
#' }
#'
#' @export
#' @md
qsf_write = function(qsf, path, format = c("qsm", "obj"), binary = TRUE)
{
  path = normalizePath(path, mustWork = FALSE)
  ext = tools::file_ext(path)

  valid_formats = c("qsm",  "qsf", "ply", "obj", "stl", "csv", "txt")
  single_file_formats = c("qsf", "obj", "ply", "stl")

  # Mode 1: Single-file export (path has a file extension)
  if (nzchar(ext))
  {
    ext_clean = tolower(ext)
    if (!ext_clean %in% valid_formats)
    {
      stop("Unsupported file extension '.", ext, "'. Must be one of: ", paste(valid_formats, collapse = ", "), call. = FALSE)
    }

    if (!ext_clean %in% single_file_formats)
    {
      stop("'.", ext, "' does not support combining several QSMs into a single file. ",
           "Single-file output is only supported for: ", paste(single_file_formats, collapse = ", "), ". ",
           "Use directory mode instead, e.g. qsf_write(qsf, \"", tools::file_path_sans_ext(basename(path)), "\", format = \"", ext_clean, "\").",
           call. = FALSE)
    }

    qsf_write_cpp(qsf, path, ext_clean, binary)
    return(invisible(TRUE))
  }

  # Mode 2: Directory export (one file per QSM inside path/format/)
  if (!is.character(format) || length(format) == 0)
  {
    stop("'format' must be a non-empty character vector.", call. = FALSE)
  }

  invalid_formats = setdiff(format, valid_formats)
  if (length(invalid_formats) > 0)
  {
    stop("Invalid format(s) requested: ", paste(invalid_formats, collapse = ", "), ". Allowed formats are: ", paste(valid_formats, collapse = ", "), call. = FALSE)
  }

  for (fmt in format)
  {
    qsf_write_cpp(qsf, path, fmt, binary)
  }

  return(invisible(TRUE))
}

#' Read a QSF Manifest File
#'
#' Reads a \code{.qsf} manifest file produced by \link{qsf_write} (whenever
#' \code{"qsm"} is among the requested \code{format}s) and reconstructs the
#' Quantitative Structural Forest (QSF) it describes.
#'
#' A \code{.qsf} file is a lightweight, standalone manifest that indexes a
#' collection of \code{.qsm} files by relative path - conceptually similar to
#' a virtual raster mosaic (VRT). Currently only "virtual" manifests are
#' supported: the manifest itself holds no tree data, only relative
#' references to the \code{.qsm} files sitting next to it, which are read in
#' turn.
#'
#' @param path A string giving the path to the \code{.qsf} manifest file.
#'
#' @return A \code{qsf} object, i.e. a list of \code{qsm} objects.
#'
#' @export
#' @seealso \link{qsf_write}
#' @md
#'
#' @examples
#' \dontrun{
#' qsf_write(qsf, "forest", format = "qsm")
#' forest <- qsf_read("forest/forest.qsf")
#' }
qsf_read = function(path)
{
  path = normalizePath(path, mustWork = TRUE)
  res = qsf_read_cpp(path)
  for (i in seq_along(res)) res[[i]] <- suppressWarnings(qsm_finalize(res[[i]]))
  res = res[order(as.numeric(names(res)))]
  res = as_qsf(res)
  res
}

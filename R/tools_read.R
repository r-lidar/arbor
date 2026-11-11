# @file tools_read.R
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

#' Read a TLS/MLS point cloud
#'
#' Reads .las, .laz or .ply files into a \code{LAS} object flagged as TLS, like \code{lidR::readTLS()}, and
#' sorts the points spatially along a Morton (Z-order) curve. With \code{rlas >= 1.9.6} LAZ files are
#' decompressed in parallel. It is the recommended way to read a point cloud before processing it with arbor.
#'
#' @param files Character. Path(s) to one or several files. Several files are merged into a single point cloud.
#' @param select Character. Attributes to read. See \code{lidR::readLAS()}.
#' @param filter Character. Streaming filters applied while reading. See \code{lidR::readLAS()}.
#'
#' @return A \code{LAS} object.
#' @export
#' @examples
#' f <- system.file("extdata", "9x9.laz", package = "arbor")
#' las <- read_tls(f, select = "xyz")
read_tls = function(files, select = "*", filter = "")
{
  stopifnot(is.character(files), length(files) >= 1L)
  stopifnot(is.character(select), length(select) == 1L)
  stopifnot(is.character(filter), length(filter) == 1L)

  if (!all(file.exists(files))) stop("File(s) not found: ", paste(files[!file.exists(files)], collapse = ", "), call. = FALSE)
  ext = tolower(tools::file_ext(files))
  if (any(!ext %in% c("las", "laz", "ply"))) stop("File(s) are not las, laz or ply", call. = FALSE)

  files = normalizePath(files)

  header = rlas::read.lasheader(files[1])

  if (!rlas_has_parallel())
  {
    message_parallel_read()
    data = rlas::read.las(files, select = select, filter = filter)
  }
  else
  {
    # rlas cannot decode spatial queries (-inside...) in parallel: it would read sequentially with a warning
    parallel = !grepl("-inside", filter, fixed = TRUE)
    data = rlas::read.las(files, select = select, filter = filter, parallel = parallel)
  }

  data = fix_arbor_types(data)

  # The bounding box and the point count of the header no longer match the data
  if (nrow(data) > 0 && (length(files) > 1 || filter != ""))
    header = rlas::header_update(header, data)

  header = drop_unread_extrabytes(header, names(data))

  las = lidR::LAS(data, header, check = TRUE, index = lidR:::LIDRTLSINDEX, no_attr_name_check = TRUE)
  spatial_sort(las)
  las
}

# The arbor attributes stored as extra bytes may have been written with another type (e.g. foliage as double).
# The C++ point cloud requires the types of PointCloudDataFrame::visit_attributes(): converts them.
ARBOR_ATTRIBUTE_TYPES = c(hag = "double", pwood = "double", foliage = "integer", passage = "integer", treeID = "integer")

fix_arbor_types = function(data)
{
  for (name in intersect(names(ARBOR_ATTRIBUTE_TYPES), names(data)))
  {
    type = ARBOR_ATTRIBUTE_TYPES[[name]]
    if (typeof(data[[name]]) != type)
      data.table::set(data, j = name, value = methods::as(data[[name]], type))
  }
  data
}

# Sorts las in place along a Morton curve (C++). Spatially sorted points make the kd-tree and the knn faster.
spatial_sort = function(las)
{
  if (lidR::npoints(las) < 2) return(invisible(las))
  C_spatial_sort(las@data)
  invisible(las)
}

rlas_has_parallel = function()
{
  utils::packageVersion("rlas") >= "1.9.6"
}

# Once per session
message_parallel_read = function()
{
  if (isTRUE(getOption("arbor.parallel_read_message_shown"))) return(invisible())
  message("Parallel reading of LAZ files is available with rlas >= 1.9.6 (installed: ", utils::packageVersion("rlas"), "). Update rlas to read faster.")
  options(arbor.parallel_read_message_shown = TRUE)
  invisible()
}

# Removes from the header the extra bytes descriptions of the attributes that were not read
drop_unread_extrabytes = function(header, read)
{
  eb = header[["Variable Length Records"]][["Extra_Bytes"]][["Extra Bytes Description"]]
  if (is.null(eb)) return(header)

  eb = eb[names(eb) %in% read]
  if (length(eb) == 0)
    header[["Variable Length Records"]][["Extra_Bytes"]] = NULL
  else
    header[["Variable Length Records"]][["Extra_Bytes"]][["Extra Bytes Description"]] = eb

  header
}

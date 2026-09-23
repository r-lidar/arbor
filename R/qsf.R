# @file qsf.R
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

#' Quantitative Structural Forest
#'
#' Batch processing of QSM models with parallel execution. Trees that are flagged with
#' \link{flag_buffer} or \link{flag_small_trees} are excluded from the computation.
#' See also the [Arbor book](https://r-lidar.github.io/arbor_book/) for more details.
#'
#' @param las A point cloud with semantic and instance segmentation computed.
#' @param min_height numeric. Default: any instance higher than 2 m generates a QSM. This is
#' pretty equivalent to \link{flag_small_trees}, but forces a height limit even if the user
#' did not call \link{flag_small_trees}, thus ensuring that no QSM lower than 2 m is computed.
#' @param params list. See \link{parameters}.
#' @return A qsf object.
#'
#' @export
#' @md
#' @seealso \link{qsm}
qsf <- function(las, min_height = 2, params = arbor_parameters_default)
{
  if (!"UserData" %in% names(las)) las@data$UserData <- ARBORTREE
  res <- qsf_cpp(las@data, min_height, params)
  for (i in seq_along(res)) res[[i]] <- suppressWarnings(qsm_finalize(res[[i]]))
  res <- res[order(as.numeric(names(res)))]
  res <- as_qsf(res)
  res
}

as_qsf <- function(x)
{
  if (!is.list(x)) {
    stop("`x` must be a list.", call. = FALSE)
  }

  if (length(x) > 0 && !all(vapply(x, inherits, logical(1), what = "qsm"))) {
    stop("All elements of `x` must be QSM objects (class 'qsm').", call. = FALSE)
  }

  class(x) <- c("qsf", class(x))
  x
}

#' Subset a qsf object
#'
#' Subsets a `qsf` object while preserving its class.
#'
#' @param x A `qsf` object.
#' @param i Index specifying the elements to extract.
#' @param ... Additional arguments passed to `[`.
#'
#' @return A `qsf` object containing the selected QSMs.
#'
#' @export
`[.qsf` <- function(x, i, ...)
{
y <- NextMethod("[")
class(y) <- class(x)
y
}

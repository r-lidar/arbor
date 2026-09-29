# @file crs.R
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

#' Get or set the projection of a QSM and QSF objects
#'
#' Get or set the projection of a QSM and QSF objects
#'
#' @param x a QSM or QSF
#' @param value see sf::st_crs
#' @param ... Unused.
#'
#' @export
#' @importFrom sf st_crs
#' @importFrom sf st_crs<-
#' @name st_crs
#' @md
#' @examples
#' f <- system.file("extdata", "qsm/tree_12.qsm", package="arbor")
#' qsm <- qsm_read(f)
#' sf::st_crs(qsm) <- 32619
#'
#' f <- system.file("extdata", "oak-plantation.qsf", package="arbor")
#' oak <- qsf_read(f)
#' sf::st_crs(oak) <- 32619
NULL

#' @export
#' @rdname st_crs
st_crs.qsm = function(x, ...)
{
  crs <- attr(x, "crs")
  if (is.null(crs)) crs <- sf::NA_crs_
  else if (crs == "") crs <- sf::NA_crs_
  else crs <- sf::st_crs(crs)
  crs
}

#' @export
#' @rdname st_crs
st_crs.qsf = function(x, ...)
{
  sf::st_crs(x[[1]])
}


#' @export
#' @rdname st_crs
`st_crs<-.qsm` = function(x, value) { attr(x, "crs") = clean_crs(sf::st_crs(value)) ; as_qsm(x) }

#' @export
#' @rdname st_crs
`st_crs<-.qsf` = function(x, value)
{
  wkt = clean_crs(sf::st_crs(value))
  y <- lapply(x, function(y) { attr(y, "crs") = wkt ; y })
  as_qsf(y)
}

clean_crs <- function(crs)
{
  if (inherits(crs, "crs")) crs <- crs$wkt
  if (is.na(crs)) return("")
  crs <- gsub("\\s+", " ", crs)
  crs <- gsub("\\s*([\\[\\],])\\s*", "\\1", crs)
  trimws(crs)
}


#' Transform the coordinates of a QSM or a QSF to a new CRS
#'
#' Method for [sf::st_transform()] applied to `qsm` and `qsf` objects. The start
#' and end points of every cylinder are transformed to the target coordinate
#' reference system. The CRS of the returned object is updated accordingly.
#'
#' @param x An object of class `qsm` or `qsf`.
#' @param crs Target coordinate reference system: an object of class `crs`, or
#'   anything accepted by [sf::st_crs()]
#' @param ... Additional arguments passed to [sf::st_transform()].
#'
#' @return An object of the same class as `x` with the columns transformed, and with
#'   its CRS set to `crs`.
#'
#' @examples
#' f <- system.file("extdata", "qsm/tree_12.qsm", package="arbor")
#' qsm <- qsm_read(f)
#' sf::st_crs(qsm) <- 32619
#' qsm2 <- sf::st_transform(qsm, 2949)
#' sf::st_crs(qsm2)
#' qsm2
#'
#' f <- system.file("extdata", "oak-plantation.qsf", package="arbor")
#' oak <- qsf_read(f)
#' sf::st_crs(oak) <- 32619
#' oak2 <- sf::st_transform(oak, 2949)
#' sf::st_crs(oak2)
#' @seealso [sf::st_transform()], [sf::st_crs()]
#' @name st_transform
#' @importFrom sf st_transform
NULL

#' @rdname st_transform
#' @method st_transform qsm
#' @export
st_transform.qsm <- function(x, crs, ...)
{
  src_crs <- sf::st_crs(x)

  if (is.na(src_crs))
    stop("Cannot transform a qsm with a missing CRS. Set it first with `sf::st_crs(x) <- `.", call. = FALSE)

  if (missing(crs))
    stop("Argument 'crs' is missing.", call. = FALSE)

  dst_crs <- sf::st_crs(crs)

  if (is.na(dst_crs))
    stop("The target 'crs' is missing or invalid.", call. = FALSE)

  if (isTRUE(sf::st_is_longlat(dst_crs)))
    warning("Transforming to a geographic CRS: radius and length attributes are not converted and remain in the original units.", call. = FALSE)

  n <- nrow(x)

  if (n > 0L)
  {
    # Start and end points are stacked so that a single transformation is done
    pts <- data.frame(
      X = c(x[["startX"]], x[["endX"]]),
      Y = c(x[["startY"]], x[["endY"]]),
      Z = c(x[["startZ"]], x[["endZ"]]))

    pts <- sf::st_as_sf(pts, coords = c("X", "Y", "Z"), crs = src_crs)
    pts <- sf::st_transform(pts, dst_crs, ...)
    xyz <- sf::st_coordinates(pts)
    xyz <- round(xyz, 3)

    i <- seq_len(n)
    x[["startX"]] <- xyz[i, "X"]
    x[["startY"]] <- xyz[i, "Y"]
    x[["startZ"]] <- xyz[i, "Z"]
    x[["endX"]]   <- xyz[n + i, "X"]
    x[["endY"]]   <- xyz[n + i, "Y"]
    x[["endZ"]]   <- xyz[n + i, "Z"]
  }

  sf::st_crs(x) <- dst_crs
  return(x)
}

#' @rdname st_transform
#' @method st_transform qsf
#' @export
st_transform.qsf <- function(x, crs, ...)
{
  if (missing(crs))
    stop("Argument 'crs' is missing.", call. = FALSE)

  # x[] <- keeps the attributes and the class of the list
  x[] <- lapply(x, st_transform.qsm, crs = crs, ...)
  return(x)
}

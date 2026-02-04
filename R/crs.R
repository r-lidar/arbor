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
NULL

#' @export
#' @rdname st_crs
st_crs.qsm = function(x, ...) { attr(x, "crs") }

#' @export
#' @rdname st_crs
st_crs.qsf = function(x, ...) { attr(x, "crs") }


#' @export
#' @rdname st_crs
`st_crs<-.qsm` = function(x, value) { attr(x, "crs") = value ; x }

#' @export
#' @rdname st_crs
`st_crs<-.qsf` = function(x, value) { attr(x, "crs") = value ; x }

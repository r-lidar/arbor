#' Subset a QSM or QSF
#'
#' Subset a QSM or QSF retaining only the stem, the merchantable or other options.\cr\cr
#' **`qsm_stem()`** retains only the stem of the QSM (branching_order = 1). Applying the function
#' to a QSF is like a loop on each QSM.\cr\cr
#' **`qsm_merchantable()`** retains only the merchantable parts of the QSM (radius > merchantable_radius
#' with additional graph validity constrains). Applying the function to a QSF is like a loop on each QSM.\cr\cr
#' **`qsf_merchantable()`** retains only the merchantable trees of the QSF (DBH > merchantable_radius). It
#' is not equivalent to `qsm_merchantable()` as it does not subset the QSMs, it subsets the QSF.
#' The entire trees are retains.
#'
#' @param qs A QSM or QSF
#' @param qsm A QSM
#' @param qsf A QSF
#'
#' @name qsm_subset
#' @rdname qsm_subset
#' @md
NULL

#' @export
#' @rdname qsm_subset
#' @export
#' @exportMethod qsm_stem
qsm_stem = function(qs)
{
  UseMethod("qsm_stem")
}

#' @export
qsm_stem.qsm = function(qs)
{
  qsm <- qsm_stem_cpp(qs)
  as_qsm(qsm)
}

#' @export
qsm_stem.qsf = function(qs)
{
  lapply(qs, qsm_stem) |> as_qsf()
}

#' @export
#' @param merchantable_radius The radius considered merchantable. Default 4.5 cm (9 cm diameter) which
#' is the standard in Canada.
#' @rdname qsm_subset
#' @exportMethod qsm_merchantable
qsm_merchantable = function(qs, merchantable_radius = 0.045)
{
  UseMethod("qsm_merchantable")
}

#' @export
qsm_merchantable.qsm = function(qs, merchantable_radius = 0.045)
{
  qsm = qsm_merchantable_cpp(qs, merchantable_radius)
  as_qsm(qsm)
}

#' @export
qsm_merchantable.qsf = function(qs, merchantable_radius = 0.045)
{
  qsf <- lapply(qs, qsm_merchantable)
  qsf <- Filter(function(x) nrow(x) > 0, qsf)
  as_qsf(qsf)
}

#' @export
#' @rdname qsm_subset
qsf_merchantable = function(qsf, merchantable_radius = 0.045)
{
  keep <- vapply(qsf, function(x)
  {
      dbh <- qsm_dbh(x)$dbh
      dbh / 2 > merchantable_radius
  }, logical(1))

  qsf[keep] |> as_qsf()
}

#' @export
#' @param stump_height The height to remove from the root. Default 25 cm
#' @rdname qsm_subset
qsm_nostump = function(qsm, stump_height = 0.25)
{
  #qsm = qsm_nostump_cpp(qsm, stump_height)
  #as_qsm(qsm)
}


#' Quantitative Structural Forest
#'
#' Batch processing of QSM models with parallel execution.
#'
#' @param las A point cloud with semantic and instance segmentation computed
#' @param params list See \link{parameters}.
#' @return A qsf object
#'
#' @export
#' @seealso  \link{qsm}
qsf <- function(las, params = default_arbor_parameters)
{
  params <- evaluate_penalty(params)
  res = qsf_cpp(las@data, params)
  for(i in seq_along(res)) res[[i]] <- qsm_finalize(res[[i]])
  res = set_qsf_class(res)
  res
}

set_qsf_class <- function(x)
{
  class(x) <- c("qsf", class(x))
  x
}

#' QSF log
#'
#' Use qsf_log() after qsf() to get the logs
#'
#' @param qsf qsf
#' @export
qsf_log = function(qsf)
{
  attr(qsf, "log")
}

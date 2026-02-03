#' @method plot qsf
#' @export
#' @rdname plot
plot.qsf = function(qsf, ...)
{
  qsf <- data.table::rbindlist(qsf)
  plot_qsm(qsf, ...)
}


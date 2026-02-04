#' @method plot qsf
#' @export
#' @rdname plot
plot.qsf = function(x, ...)
{
  qsf <- data.table::rbindlist(x)
  plot_qsm(qsf, ...)
}


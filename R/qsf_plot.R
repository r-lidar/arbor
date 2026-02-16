#' @method plot qsf
#' @export
#' @rdname plot
plot.qsf = function(x, ...)
{
  qsf <- data.table::rbindlist(x)
  qsf$cyl_ID = 1:nrow(qsf)
  plot_qsm(qsf, ...)
}


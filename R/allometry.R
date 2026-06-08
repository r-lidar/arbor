#' Available allometric models
#'
#' Show available allometric models
#'
#' @return
#' A data.frame containing available models and their reference URLs:
#' \describe{
#'   \item{model}{Model name}
#'   \item{url}{DOI or reference link}
#' }
#'
#' @examples
#' available_allometries()
#'
#' @export
available_allometries <- function()
{
  models <- data.frame(
    model = c(
      "CostaCysneiros2020",
      "Griese2025",
      "Chenge2020"
    ),
    country = c(
      "Brazil",
      "Germany",
      "Nigeria"
    ),
    color = c(
      "#1b9e77",
      "#d95f02",
      "#7570b3"
    ),
    url = c(
      "https://doi.org/10.1139/cjfr-2020-0060",
      "https://doi.org/10.1038/s41597-025-06421-7",
      "https://doi.org/10.1016/j.tfp.2020.100051"
    ),
    stringsAsFactors = FALSE
  )

  filters <- setNames(
    replicate(nrow(models),
              function(x) { x$DBH <- x$DBH * 100; x },
              simplify = FALSE),
    models$model
  )

  curves <- Map(
    function(name, fun) fun(arbor:::allometry(name)),
    models$model,
    filters
  )

  plot(
    curves[[1]],
    col  = models$color[1],
    lwd  = 2,
    xlab = "DBH (cm)",
    ylab = "Height",
    type = "l",
    xlim = c(0, 50)
  )

  for (i in seq_along(curves)[-1]) {
    lines(
      curves[[i]],
      col = models$color[i],
      lwd = 2
    )
  }

  legend(
    "bottomright",
    legend = sprintf("%s (%s)", models$model, models$country),
    col = models$color,
    lwd = 2,
    bty = "n"
  )

  abline(v = 6, lty = 3)

  library_table <- models[, c("model", "url")]

  return(library_table)
}

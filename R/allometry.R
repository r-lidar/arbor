#' Available allometric models
#'
#' Show available allometric models
#'
#' @return
#' A data.frame containing available models and their references:
#' \describe{
#'   \item{model}{Model name}
#'   \item{authors}{Authors of the model}
#'   \item{country}{Country where the model was calibrated}
#'   \item{reference}{DOI or reference link}
#' }
#'
#' @examples
#' available_allometries()
#'
#' @export
available_allometries <- function()
{
  models <- allometry_database_cpp()

  palette <- c("#1b9e77", "#d95f02", "#7570b3", "#8e129e", "#e7298a", "#66a61e", "#e6ab02", "#a6761d")
  colors  <- rep_len(palette, nrow(models))

  curves <- lapply(models$model, function(name) {
    x <- allometry(name)
    x$DBH <- x$DBH * 100
    x
  })

  plot(
    curves[[1]],
    col  = colors[1],
    lwd  = 2,
    xlab = "DBH (cm)",
    ylab = "Height",
    type = "l",
    xlim = c(0, 50)
  )

  for (i in seq_along(curves)[-1]) {
    graphics::lines(
      curves[[i]],
      col = colors[i],
      lwd = 2
    )
  }

  graphics::legend(
    "bottomright",
    legend = sprintf("%s (%s)", models$model, models$country),
    col = colors,
    lwd = 2,
    bty = "n"
  )

  return(models)
}

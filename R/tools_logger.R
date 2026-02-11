logger <- function(
    msg,
    level = NULL,     # NULL = no level printed
    time = TRUE,
    file = NULL,
    ...
) {
  # ANSI color codes (base R only)
  color <- "\033[0m"
  if (!is.null(level))
  {
    color <- switch(
      level,
      "INFO"  = "\033[34m",  # blue
      "WARN"  = "\033[33m",  # yellow
      "ERROR" = "\033[31m",  # red
      "DEBUG" = "\033[36m",  # cyan
      ""
    )
  }
  reset <- "\033[0m"

  timestamp <- if (time) sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")) else ""
  level_tag <- if (!is.null(level)) sprintf("[%s] ", level) else ""

  message <- paste(msg, ..., sep = " ")
  out <- paste0(timestamp, level_tag, message)

  if (is.null(file)) {
    cat(color, out, reset, "\n", sep = "")
  } else {
    cat(out, "\n", file = file, append = TRUE)
  }

  invisible(out)
}

log_time <- function(label, expr) {
  t0 <- Sys.time()
  logger(paste("START", label))
  on.exit({
    dt <- difftime(Sys.time(), t0, units = "secs")
    logger(sprintf("END %s (%.2f s)", label, as.numeric(dt)))
  })
  force(expr)
}

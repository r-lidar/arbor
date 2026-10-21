mem_monitor <- function(output = file.path(tempdir(), "mem_trace.csv"), pulse = 1) {
  script <- system.file("bash", "mem_monitor.sh", package = "arbor")
  if (script == "") stop("mem_monitor.sh not found in package 'arbor'")

  phase_file <- tempfile("mem_phase_", fileext = ".txt")
  writeLines("start", phase_file)

  system2("setsid",
          c("bash", shQuote(script), Sys.getpid(), shQuote(output),
            as.character(pulse), shQuote(phase_file)),
          wait = FALSE, stdout = FALSE, stderr = FALSE)

  invisible(phase_file)
}

phase <- function(x, file) writeLines(x, file)

mem_monitor_stop <- function(file) invisible(unlink(file))

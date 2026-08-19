#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lidR)
  library(arbor)
  library(terra)
  library(tools)
})

# ----------
# Report
# ----------

run_report = function(args)
{
  usage_report = function()
  {
    cat("
Usage:
  arbor report <input_dir> <output_pdf>
")
    quit(save = "no", status = 0)
  }

  if (length(args) < 2 || arbor:::has_flag(args, "-h") || arbor:::has_flag(args, "--help")) {
    usage_report()
  }

  input_dir <- normalizePath(args[1])
  output_pdf <- normalizePath(args[2], mustWork = FALSE)

  if (!dir.exists(input_dir)) {
    stop("Input directory does not exist: ", input_dir)
  }

  rmarkdown::render(
    input  = system.file("bash", "report_template.Rmd", package="arbor"),
    output_file = basename(output_pdf), # only filename
    output_dir = dirname(output_pdf), # directory to write to
    params = list(
      idir = input_dir
    ),
    quiet = FALSE
  )
}


# ------------------------------------------------------------
# Main Dispatcher
# ------------------------------------------------------------

all_args <- commandArgs(trailingOnly = TRUE)

if (length(all_args) == 0) {
  cmd_usage()
}

command  <- all_args[1]
sub_args <- all_args[-1]

t0 = arbor:::tic()

if (command == "segment") {
  arbor:::cmd_segment(sub_args)
} else if (command == "qsf") {
  arbor:::cmd_qsf(sub_args)
} else if (command == "report") {
  run_report(sub_args)
} else if (command == "blender") {
  arbor:::cmd_blender(sub_args)
} else if (command %in% c("-h", "--help")) {
  arbor:::cmd_usage()
} else {
  cat(paste0("Unknown command: '", command, "'\n"))
  arbor:::cmd_usage()
}

arbor:::toc(t0, space = "")

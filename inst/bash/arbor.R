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
  arbor report <input.qsf>
    
Mandatory:
  <input.qsf>             Input QSF file
    
Options:
  -epsg 1234              EPSG code to georeference the dataset if missing
")
    quit(save = "no", status = 0)
  }

  if (length(args) < 1 || arbor:::has_flag(args, "-h") || arbor:::has_flag(args, "--help")) {
    usage_report()
  }

  epsg <- as.numeric(arbor:::get_arg(args, "-epsg", 0))
  input_qsf <- normalizePath(args[1])

  if (!file.exists(input_qsf)) {
    stop("Input directory does not exist: ", input_qsf)
  }

  rmarkdown::render(
    input  = system.file("bash", "report_template_html.Rmd", package="arbor"),
    output_file = paste0(tools::file_path_sans_ext(input_qsf), ".html"),
    output_dir = dirname(input_qsf),   # directory to write to
    params = list(
      input = input_qsf,
      epsg = epsg
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

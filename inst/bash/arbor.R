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

# Options known by `arbor report`
REPORT_VALUE_FLAGS <- c("-o", "-target_epsg")
REPORT_BOOL_FLAGS  <- c("-h", "--help")

run_report = function(args, strict = TRUE)
{
  usage_report = function()
  {
    cat("
Usage:
  arbor report <input.qsf>

Mandatory:
  <input.qsf>             Input QSF file

Options:
  -o                      Path to exported html file
  -target_epsg 1234       EPSG code to georeference the dataset if missing
")
    quit(save = "no", status = 0)
  }

  if (length(args) < 1 || arbor:::has_flag(args, "-h") || arbor:::has_flag(args, "--help")) {
    usage_report()
  }

  if (strict) arbor:::check_args(args, REPORT_VALUE_FLAGS, REPORT_BOOL_FLAGS)

  target_epsg <- as.numeric(arbor:::get_arg(args, "-target_epsg", 0))
  output      <- arbor:::get_arg(args, "-o", "")
  input_qsf   <- normalizePath(args[1])

  if (output == "")
  {
    if (dir.exists(input_qsf)) {
      # Places <folder_name>.html INSIDE the directory
      output <- file.path(input_qsf, paste0(basename(input_qsf), ".html"))
    } else {
      # Standard behavior for file paths
      output <- paste0(tools::file_path_sans_ext(input_qsf), ".html")
    }
  }

  output = normalizePath(output, mustWork = FALSE)

  if (!file.exists(input_qsf)) {
    stop("Input directory does not exist: ", input_qsf)
  }

  rmarkdown::render(
    input  = system.file("bash", "report_template_html.Rmd", package="arbor"),
    output_file = basename(output),
    output_dir = dirname(output),   # directory to write to
    params = list(
      input = input_qsf,
      target_epsg = target_epsg,
      output = output
    ),
    quiet = FALSE
  )
}

# ----------
# Pipeline (segment + qsf + report)
# ----------

run_pipeline = function(args)
{
  usage_pipeline = function()
  {
    cat("
Usage:
  arbor pipeline <input.las> [options]

Runs segment -> qsf -> report as one pass. Any option accepted by
'arbor segment' or 'arbor qsf' can be passed here as-is and is forwarded
to that step -- run 'arbor segment -h' / 'arbor qsf -h' for the full list.
-target_epsg is forwarded to the report. -o is only used by the qsf step; the
report is written next to the .qsf file.

Mandatory:
  <input.las>              Input LAS/LAZ file (or a trees file if -skip-segment is used)

Pipeline control:
  -skip-segment            <input.las> is already a segmented trees file; start at the qsf step
  -skip-report             Stop after the qsf step, do not render the report
  -h, --help               Show this help
")
    quit(save = "no", status = 0)
  }

  if (length(args) < 1 || arbor:::has_flag(args, "-h") || arbor:::has_flag(args, "--help")) {
    usage_pipeline()
  }

  input <- args[1]
  if (startsWith(input, "-")) {
    stop("First argument must be the input file.")
  }
  if (!file.exists(input)) {
    stop("Input file does not exist: ", input)
  }
  input <- normalizePath(input, mustWork = TRUE)

  # Everything after the input file is forwarded to the steps. Each step only
  # reads the flags it knows about, so the whole command line is validated once
  # here against the union of the flags of all steps (each step then runs with
  # strict = FALSE). The flag lists are owned by each command.
  rest_args <- args[-1]

  arbor:::check_args(args,
    value_flags = c(arbor:::SEGMENT_VALUE_FLAGS, arbor:::QSF_VALUE_FLAGS, REPORT_VALUE_FLAGS),
    bool_flags  = c(arbor:::SEGMENT_BOOL_FLAGS, arbor:::QSF_BOOL_FLAGS, REPORT_BOOL_FLAGS, "-skip-segment", "-skip-report"))

  skip_segment <- arbor:::has_flag(rest_args, "-skip-segment")
  skip_report  <- arbor:::has_flag(rest_args, "-skip-report")

  # ---- 1. Segment ----
  if (skip_segment)
  {
    cat("\n=== [Pipeline 1/3] Segmentation skipped (-skip-segment) ===\n")
    seg_res <- list(dir = dirname(input), trees = input)
  }
  else
  {
    cat("\n=== [Pipeline 1/3] Segmentation ===\n")
    seg_res <- arbor:::cmd_segment(c(input, rest_args), strict = FALSE)
    gc()
  }

  if (is.null(seg_res) || is.na(seg_res$trees))
  {
    stop("Pipeline stopped: no trees file available to feed into the QSF step")
  }

  # ---- 2. QSM ----
  qsf_args <- c(seg_res$trees, rest_args)

  cat("\n=== [Pipeline 2/3] QSM ===\n")
  qsf_res <- arbor:::cmd_qsf(qsf_args, strict = FALSE)
  gc()

  if (is.null(qsf_res) || is.na(qsf_res$qsf_file))
  {
    stop("Pipeline stopped: the QSF step did not produce a .qsf file for the report step.")
  }

  # ---- 3. Report ----
  if (skip_report)
  {
    cat("\n=== [Pipeline 3/3] Report skipped (-skip-report) ===\n")
    return(invisible(list(segment = seg_res, qsf = qsf_res, report = NA_character_)))
  }

  # '-o' is a directory for the qsf step but an html file for the report: do not forward it
  report_args <- qsf_res$qsf_file
  target_epsg <- arbor:::get_arg(rest_args, "-target_epsg")
  if (!is.null(target_epsg)) report_args <- c(report_args, "-target_epsg", target_epsg)

  cat("\n=== [Pipeline 3/3] Report ===\n")
  report_path <- run_report(report_args, strict = FALSE)

  invisible(list(segment = seg_res, qsf = qsf_res, report = report_path))
}


# ------------------------------------------------------------
# Main Dispatcher
# ------------------------------------------------------------

all_args <- commandArgs(trailingOnly = TRUE)

if (length(all_args) == 0) {
  arbor:::cmd_usage()
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
} else if (command == "pipeline") {
  run_pipeline(sub_args)
} else if (command %in% c("-h", "--help")) {
  arbor:::cmd_usage()
} else {
  cat(paste0("Unknown command: '", command, "'\n"))
  arbor:::cmd_usage()
}

arbor:::toc(t0, space = "")

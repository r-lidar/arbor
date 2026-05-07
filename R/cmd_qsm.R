cmd_qsm <- function(args) {

  # --- QSM Usage ---
  usage_qsm <- function() {
    cat("
Usage:
  arbor qsm <input> [options]

Mandatory:
  <input>           Input file

Options:
  -o, --output <dir> Output directory [default: input directory]
  -overwrite         Overwrite existing files
  -csv               Export CSV [default: on if nothing provided]
  -obj               Export OBJ [default: on if nothing provided]
  -ply               Export PLY

Example:
  arbor qsm tree.laz -csv -obj
")
    quit(save = "no", status = 0)
  }

  if (length(args) == 0 || has_flag(args, "-h") || has_flag(args, "--help")) {
    usage_qsm()
  }

  # --- Parse Inputs ---
  # Prioritize positional argument for input
  ifiles <- NULL
  if (!startsWith(args[1], "-")) {
    ifiles <- args[1]
  } else {
    # Fallback to flags if user uses legacy style
    ifiles <- get_arg(args, "-i")
    if (is.null(ifiles)) ifiles <- get_arg(args, "-input")
  }

  if (is.null(ifiles)) fail("Missing input file or folder")
  if (!file.exists(ifiles)) fail(paste("Input does not exist:", ifiles))
  ifiles <- normalizePath(ifiles, mustWork = TRUE)

  # Output dir
  odir <- file.path(dirname(ifiles), "qsm")

  # Override output if flag present
  out_flag <- get_arg(args, "-o")
  if (is.null(out_flag)) out_flag <- get_arg(args, "--output")
  if (!is.null(out_flag)) odir <- normalizePath(out_flag, mustWork = FALSE)

  # Formats
  formats <- character()
  if (has_flag(args, "-csv")) formats <- c(formats, "csv")
  if (has_flag(args, "-obj")) formats <- c(formats, "obj")
  if (has_flag(args, "-ply")) formats <- c(formats, "ply")
  if (length(formats) == 0) formats <- c("csv", "obj") # Default

  overwrite <- has_flag(args, "-overwrite")

  cat("
================= Arbor QSM module =================
Input       :", paste(ifiles, collapse = ", "), "
Output      :", odir, "
Settings
  Formats   :", paste(formats, collapse = ", "), "
=====================================================
")

  # res <- qsf(
  #   input     = ifiles,
  #   odir      = odir,
  #   formats   = formats,
  #   overwrite = overwrite,
  #   ncores    = ncores
  # )
  #
  # log <- qsf_log(res)

  #for (i in which(has_msg))
  #{
  #  if (isTRUE(res$success[i])) {
  #    cat("\033[33m", "WARNING | ", res$name[i], " | ", res$message[i], "\033[0m\n", sep = "")
  #  } else {
  #    cat("\033[31m", "ERROR   | ", res$name[i], " | ", res$message[i], "\033[0m\n", sep = "")
  #  }
  #}
}

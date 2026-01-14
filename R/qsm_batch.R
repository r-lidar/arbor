#' Batch processing of QSM models with parallel execution
#'
#' This function processes multiple LAS files using \link{qsm} in parallel.
#' It supports exporting CSV/OBJ/PLY QSM outputs, optional overwriting of
#' existing files, and returns detailed logs for each processed file.
#'
#' @param ifiles Character vector of input file paths, or a single directory containing LAS/LAZ files.
#' @param odir Output directory where `csv/`, `obj/`, and `ply/` folders will be created.
#' @param formats Character; write QSM results in which formats? See \link{qsm_write}.
#' @param overwrite Logical; overwrite existing outputs. Default: FALSE. Avoid recomputing already
#' existing QSMs
#' @param ncores Number of CPU cores to use.
#' @param ... Additional arguments passed to \link{qsm}.
#'
#' @return A data.table of logs, one line per input file. Each log entry contains:
#'   \itemize{
#'     \item `file`: input file path
#'     \item `name`: basename without extension
#'     \item `success` TRUE/FALSE
#'     \item `status`: "success", "skipped_existing", "failed_readLAS", "failed_qsm"
#'     \item `error`: error message if any
#'   }
#'
#' @export
#' @seealso  \link{qsm}
qsm_batch = function(
    ifiles,
    odir,
    formats = c("csv", "obj"),
    overwrite = FALSE,
    ncores = parallel::detectCores() / 2,
    ...
) {
  ti <- tic()

  # 1. Validation and Directory Setup
  formats <- unique(tolower(formats))
  if (length(formats) == 0) stop("At least one format must be specified.")

  odir <- normalizePath(odir, mustWork = FALSE)

  # Create a subfolder for every requested format automatically
  format_dirs <- stats::setNames(file.path(odir, formats), formats)
  for (d in format_dirs) {
    if (!dir.exists(d))
      dir.create(d, recursive = TRUE)
  }

  # 2. Handle input files/directories
  if (length(ifiles) == 1 && dir.exists(ifiles)) {
    ifiles <- list.files(
      ifiles,
      pattern = "\\.(las|laz)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
  }
  ifiles <- normalizePath(ifiles, mustWork = TRUE)
  if (length(ifiles) == 0) stop("No valid input files found.")

  # 3. Parallel Setup
  dots <- list(...)
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl))

  # Export only the necessary dynamic variables
  parallel::clusterExport(
    cl,
    varlist = c("formats", "format_dirs", "overwrite", "dots", "qsm", "qsm_write"),
    envir = environment()
  )

  # 4. Process Files
  res <- pbapply::pblapply(ifiles, cl = cl, FUN = function(f)
  {
    name <- tools::file_path_sans_ext(basename(f))

    # Generate all output paths dynamically
    out_paths <- file.path(format_dirs, paste0(name, ".", formats))
    names(out_paths) <- formats

    log <- list(
      file = f,
      name = name,
      success = TRUE,
      status = "pending",
      error = NA_character_
    )

    # Skip logic: Check if all requested files exist
    if (!overwrite && all(file.exists(out_paths))) {
      log$status <- "skipped_existing"
      return(log)
    }

    # Execute QSM
    las <- tryCatch(lidR::readLAS(f), error = function(e) e)
    if (inherits(las, "error"))
    {
      log$success <- FALSE; log$status <- "failed_readLAS"; log$error <- las$message
      return(log)
    }

    tryCatch(
    {
      # Capture logs to keep console clean during parallel execution
      utils::capture.output({
        q <- do.call(qsm, c(list(las), dots))
      }, type = "output")

      # Loop through formats and write files automatically
      # This will work as long as qsm_write supports the extension
      for (o in out_paths) qsm_write(q, o)

      log$status <- "success"
    },
    error = function(e)
    {
      log$success <- FALSE; log$status <- "failed_qsm"; log$error <- e$message
    })

    return(log)
  })

  res = data.table::rbindlist(res)

  toc(ti)

  return(res)
}

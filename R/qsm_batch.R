#' Batch processing of QSM models with parallel execution
#'
#' This function processes multiple LAS files using \link{qsm} in parallel.
#' It supports exporting CSV/OBJ/PLY QSM outputs, optional overwriting of
#' existing files, and returns detailed logs for each processed file.
#'
#' @param ifiles Character vector of input file paths, or a single directory containing LAS/LAZ files.
#' @param odir Output directory where `csv/`, `obj/`, and `ply/` folders will be created.
#' @param csv Logical; write QSM results as CSV. Default: TRUE.
#' @param obj Logical; write QSM results as OBJ. Default: TRUE.
#' @param ply Logical; write QSM results as PLY. Default: FALSE.
#' @param overwrite Logical; overwrite existing outputs. Default: FALSE. Avoid recomputing already
#' existing QSMs
#' @param ncores Number of CPU cores to use. Default: all available cores.
#' @param ... Additional arguments passed to `qsm()`.
#'
#' @return A data.table of logs, one line per input file. Each log entry contains:
#'   \itemize{
#'     \item `file`: input file path
#'     \item `name`: basename without extension
#'     \iten `success` TRUE/FALSE
#'     \item `status`: "success", "skipped_existing", "failed_readLAS", "failed_qsm"
#'     \item `error`: error message if any
#'   }
#'
#' @export
#' @seealso  \link{qsm}
qsm_batch = function(
    ifiles,
    odir,
    csv = TRUE,
    obj = TRUE,
    ply = FALSE,
    overwrite = FALSE,
    ncores = parallel::detectCores(),
    ...
)
{
  # Normalize booleans
  csv <- isTRUE(csv)
  obj <- isTRUE(obj)
  ply <- isTRUE(ply)

  # Validate overwrite
  overwrite <- isTRUE(overwrite)

  # Short-circuit
  if (!any(csv, obj, ply))
    stop("Nothing to do: csv=FALSE, obj=FALSE, ply=FALSE.")

  # Capture dots for cluster
  dots <- list(...)

  ti <- tic()

  # Build output directories
  odir     <- normalizePath(odir, mustWork = FALSE)

  ocsvdir  <- file.path(odir, "csv")
  oobjdir  <- file.path(odir, "obj")
  oplydir  <- file.path(odir, "ply")

  dirs_to_create <- c(
    if (csv) ocsvdir else NULL,
    if (obj) oobjdir else NULL,
    if (ply) oplydir else NULL
  )
  for (d in dirs_to_create) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

  # Handle directory input
  if (length(ifiles) == 1 && dir.exists(ifiles))
  {
    dirpath <- normalizePath(ifiles, mustWork = TRUE)
    ifiles <- list.files(
      dirpath,
      pattern = "\\.(las|laz)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (length(ifiles) == 0)
      stop("Directory contains no LAS/LAZ files: ", dirpath)
  }

  # List files
  ifiles <- normalizePath(ifiles, mustWork = TRUE)

  # Prepare logs
  logs <- vector("list", length(ifiles))
  names(logs) <- basename(ifiles)

  # Parallel setup
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl))

  parallel::clusterExport(
    cl,
    varlist = c(
      "csv", "obj", "ply", "ocsvdir", "oobjdir", "oplydir",
      "overwrite", "dots"
    ),
    envir = environment()
  )

  # Process files
  res <- pbapply::pblapply(ifiles, cl = cl, FUN = function(f)
  {
    name <- tools::file_path_sans_ext(basename(f))

    outcsv <- file.path(ocsvdir, paste0(name, ".csv"))
    outobj <- file.path(oobjdir, paste0(name, ".obj"))
    outply <- file.path(oplydir, paste0(name, ".ply"))

    log <- list(
      file = f,
      name = name,
      success = TRUE,
      status = "pending",
      error = NA_character_
    )

    # Skip?
    needed <- c(
      csv && !file.exists(outcsv),
      obj && !file.exists(outobj),
      ply && !file.exists(outply)
    )

    if (!overwrite && !any(needed))
    {
      log$status <- "skipped_existing"
      return(log)
    }

    # Read LAS
    las <- tryCatch(lidR::readLAS(f), error = function(e) e)

    if (inherits(las, "error"))
    {
      log$success = FALSE
      log$status <- "failed_readLAS"
      log$error <- las$message
      return(log)
    }

    # QSM execution
    tryCatch(
    {
      qsm_logs <- capture.output({
        q <- do.call(qsm, c(list(las), dots))
      }, type = "output")

      # Write outputs
      if (csv) { qsm_write(q, outcsv) }
      if (obj) { qsm_write(q, outobj) }
      if (ply) { qsm_write(q, outply) }

      log$status <- "success"
    },
    error = function(e) {
      log$success <- FALSE
      log$status  <- "failed_qsm"
      log$error   <- e$message
    })

    return(log)
  })

  res = data.table::rbindlist(res)

  toc(ti)

  return(res)
}

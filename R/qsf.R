#' Quantitative Structural Forest
#'
#' Batch processing of QSM models with parallel execution. This function processes multiple LAS files
#' using \link{qsm} in parallel. It supports exporting CSV/OBJ/PLY QSM outputs, optional overwriting
#' of existing files, and returns detailed logs for each processed file.
#'
#' @param input Character vector of input file paths, or a single directory containing LAS/LAZ files
#' or a LAS object with semantic and instance segmentation.
#' @param odir Output directory where `csv/`, `obj/`, and `ply/` folders will be created. odir can
#' be missing but in this case nothing is exported. Only a table is returned.
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
#'     \item `message`: error message if any
#'   }
#'
#' @export
#' @seealso  \link{qsm}
qsf <- function(
    input,
    odir,
    formats = c("csv", "obj"),
    overwrite = FALSE,
    ncores = parallel::detectCores() / 2,
    ...
)
{
  ti <- tic()

  if (inherits(input, "LAS"))
  {
    cat("Input is a point cloud: exporting trees in temporary files\n")
    treeID <- NULL
    dir = tempdir()
    input = sapply(unique(input$treeID), function(i)
    {
      tree <- lidR::filter_poi(input, treeID == i)
      olas <- paste0(dir, "/tree_", i, ".las")
      lidR::writeLAS(tree, olas)
    })
  }

  # ---------------------------------------------------------------------------
  # 1. Validation and directory setup
  # ---------------------------------------------------------------------------

  formats <- unique(tolower(formats))
  if (length(formats) == 0)
    stop("At least one format must be specified.")

  if (!missing(odir))
  {
    odir <- normalizePath(odir, mustWork = FALSE)

    format_dirs <- stats::setNames(file.path(odir, formats), formats)
    for (d in format_dirs)
      if (!dir.exists(d))
        dir.create(d, recursive = TRUE)
  }
  else
  {
    format_dirs = NULL
    odir = NULL
  }

  # ---------------------------------------------------------------------------
  # 2. Handle inputs
  # ---------------------------------------------------------------------------

  if (length(input) == 1 && dir.exists(input))
  {
    input <- list.files(
      input,
      pattern = "\\.(las|laz)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
  }

  input <- normalizePath(input, mustWork = TRUE)
  if (length(input) == 0)
    stop("No valid input files found.")

  dots <- list(...)

  # ---------------------------------------------------------------------------
  # 3. Future plan management (user-friendly + safe)
  # ---------------------------------------------------------------------------

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (is.null(old_plan) || inherits(old_plan, "sequential")) {
    if (ncores > 1) {
      future::plan(
        future::multisession,
        workers = as.integer(ncores)
      )
    }
  }

  # ---------------------------------------------------------------------------
  # 4. Progress-aware parallel processing
  # ---------------------------------------------------------------------------

  if (!interactive())
  {
    options(progressr.enable = TRUE)
  }
  progressr::handlers(
    progressr::handler_progress(
    format = ":spin :current/:total [:bar] :percent in :elapsed ETA: :eta",
    width = 50,
    complete = "="
  ))

  res <- progressr::with_progress(
  {
    p <- progressr::progressor(steps = length(input))

    future.apply::future_lapply(input, function(f)
    {
      on.exit(p(), add = TRUE)

      name <- tools::file_path_sans_ext(basename(f))

      out_paths = NULL
      if (!is.null(format_dirs))
      {
        out_paths <- file.path(format_dirs, paste0(name, ".", formats))
        names(out_paths) <- formats
      }

      # Initialize log
      log <- list(
        file = f,
        name = name,
        success = TRUE,
        status = "pending",
        message = NA_character_
      )

      # Temporary variable to collect warnings for this specific file
      w_msgs <- character()

      # Helper to clean up final error message (combine warnings + errors)
      format_msg <- function(warnings, err = NULL)
      {
        all_msgs <- c(warnings, err)
        if (length(all_msgs) == 0) return(NA_character_)
        return(paste(unique(all_msgs), collapse = " | "))
      }

      if (!is.null(out_paths))
      {
        if (!overwrite && all(file.exists(out_paths)))
        {
          log$status <- "skipped_existing"
          return(log)
        }
      }

      # --- Step 1: Read LAS ---
      # We use withCallingHandlers to catch warnings without stopping
      las <- tryCatch(
        withCallingHandlers(
          lidR::readLAS(f),
          warning = function(w) {
            w_msgs <<- c(w_msgs, w$message)
            invokeRestart("muffleWarning") # Suppress printing to console
          }
        ),
        error = function(e) e
      )

      if (inherits(las, "error"))
      {
        log$success <- FALSE
        log$status <- "failed_readLAS"
        log$message  <- format_msg(w_msgs, las$message)
        return(log)
      }

      # --- Step 2: QSM & Write ---
      q <- tryCatch(
      {
        utils::capture.output(
        {
          # Wrap the main calculation in withCallingHandlers as well
          withCallingHandlers(
          {
            q <- do.call(qsm, c(list(las), dots))
          },
          warning = function(w)
          {
            w_msgs <<- c(w_msgs, w$message)
            invokeRestart("muffleWarning")
          })
        }, type = "output")

        for (o in out_paths)
          qsm_write(q, o)

        log$status <- "success"

        # If there were warnings, we still mark success, but record warnings in 'error'
        if (length(w_msgs) > 0)
        {
          log$message <- format_msg(w_msgs)
        }

        q
      },
      error = function(e)
      {
        log$success <- FALSE
        log$status  <- "failed_qsm"
        log$message   <- format_msg(w_msgs, e$message)
      })

      attr(q, "name") = name
      attr(q, "log") = log
      q
    },
    future.seed = TRUE)
  })

  log <- lapply(res, function(x) attr(x, "log"))
  log <- data.table::rbindlist(log)

  names <- sapply(res, function(x) attr(x, "name"))
  names(res) <- names
  res <- set_qsf_class(res)
  attr(res, "log") = log

  toc(ti)
  res
}

set_qsf_class <- function(x)
{
  class(x) <- c("qsf", class(x))
  x
}

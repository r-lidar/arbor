# @file qsf.R
# Project: Arbor
#
# Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' Quantitative Structural Forest
#'
#' Batch processing of QSM models with parallel execution.
#'
#' @param las A point cloud with semantic and instance segmentation computed
#' @param min_height numeric. Default: any instance higher than 2 m generates a QSMs
#' @param params list See \link{parameters}.
#' @return A qsf object
#'
#' @export
#' @seealso  \link{qsm}
qsf <- function(las, min_height = 2, params = arbor_parameters_default)
{
  if (!"UserData" %in% names(las)) las@data$UserData <- ARBORTREE
  res <- qsf_cpp(las@data, min_height, params)
  for (i in seq_along(res)) res[[i]] <- suppressWarnings(qsm_finalize(res[[i]]))
  res <- res[order(as.numeric(names(res)))]
  res <- as_qsf(res)
  res
}

as_qsf <- function(x)
{
  if (!is.list(x)) {
    stop("`x` must be a list.", call. = FALSE)
  }

  if (length(x) > 0 && !all(vapply(x, inherits, logical(1), what = "qsm"))) {
    stop("All elements of `x` must be QSM objects (class 'qsm').", call. = FALSE)
  }

  class(x) <- c("qsf", class(x))
  x
}

#' Subset a qsf object
#'
#' Subsets a `qsf` object while preserving its class.
#'
#' @param x A `qsf` object.
#' @param i Index specifying the elements to extract.
#' @param ... Additional arguments passed to `[`.
#'
#' @return A `qsf` object containing the selected QSMs.
#'
#' @export
`[.qsf` <- function(x, i, ...)
{
y <- NextMethod("[")
class(y) <- class(x)
y
}


#' Registry of known QSF log codes
#'
#' Internal registry mapping a log `code` (e.g. `"W0"`) to its default short
#' label. The leading letter of a code identifies its `type`: `W` = warning,
#' `M` = message, `E` = error. This registry is used as a fallback to label
#' coded messages that do not repeat their label in the message text, and to
#' document the codes produced by the underlying QSM reconstruction.
#'
#' @keywords internal
.arbor_qsf_codes <- c(
  W0 = "No wood point",
  W1 = "QSM Generation Failed",
  W2 = "No valid measure",
  W3 = "Small tree allometry",
  W4 = "Broken tree",
  W5 = "Diameter anomaly"
)

#' Get the type (warning/message/error) associated with a log code
#' @param code character. A code such as `"W0"`.
#' @keywords internal
.qsf_code_type = function(code)
{
  if (is.na(code)) return(NA_character_)

  letter <- substr(code, 1, 1)
  switch(letter,
    W = "warning",
    M = "message",
    E = "error",
    NA_character_
  )
}

#' Parse a single QSM log message
#'
#' Parses a message of the form `"[W0] [Short label] Long explanation"` (coded,
#' preferred) or the legacy `"[Short label] Long explanation"` (uncoded). This
#' is the single source of truth for message parsing, used by both
#' \link{qsf_log} and \link{qsm_message}.
#'
#' @param msg character(1). A single log message.
#' @return A list with `code`, `type`, `label` and `message`.
#' @keywords internal
.qsf_parse_message = function(msg)
{
  if (length(msg) == 0 || is.na(msg))
    return(list(code = NA_character_, type = NA_character_, label = NA_character_, message = msg))

  code_pattern  <- "^\\s*\\[([EWM][0-9]+)\\]\\s*(.*)$"
  label_pattern <- "^\\s*\\[([^]]+)\\]\\s*(.*)$"

  code <- NA_character_
  rest <- msg

  if (grepl(code_pattern, msg))
  {
    code <- sub(code_pattern, "\\1", msg)
    rest <- sub(code_pattern, "\\2", msg)
  }

  label <- NA_character_
  if (grepl(label_pattern, rest))
  {
    label <- sub(label_pattern, "\\1", rest)
  }
  else if (!is.na(code) && code %in% names(.arbor_qsf_codes))
  {
    label <- unname(.arbor_qsf_codes[code])
  }

  list(code = code, type = .qsf_code_type(code), label = label, message = msg)
}

#' QSF log
#'
#' Use qsf_log after \link{qsf} to get a structured summary of the log messages
#' recorded during the QSM reconstruction of every tree in a `qsf` object.
#'
#' @param qsf qsf
#'
#' @return
#' A named list, one element per distinct log entry. Elements are named by
#' their `code` (e.g. `"W0"`) when the underlying message follows the code
#' convention (see Codes section below), or fall back to the short label text
#' for legacy, uncoded messages, preserving the historical behavior. Each
#' element is a list with:
#' \itemize{
#'   \item{`code`: the log code (e.g. `"W0"`), or `NA` for legacy, uncoded messages.}
#'   \item{`type`: `"warning"`, `"message"`, `"error"`, or `NA` if unknown.}
#'   \item{`label`: the short, human readable label extracted from the message.}
#'   \item{`message`: an example of the full message text.}
#'   \item{`index`: the indices, within `qsf`, of the trees carrying this log entry.}
#'   \item{`treeID`: the tree IDs of the trees carrying this log entry.}
#' }
#'
#' @section Codes:
#' Messages recorded on a QSM may start with a code `[<L><n>]`, where `<L>` is
#' one of `W` (warning), `M` (message) or `E` (error) and `<n>` is a numeric
#' identifier, immediately followed by the historical `[Short label]` and the
#' free text explanation, e.g. `"[W0] [No wood point] This tree has no point
#' labelled as wood"`. This is fully backward compatible: messages recorded
#' before this convention was introduced have no code (`code = NA`) and are
#' grouped by their short label instead, exactly as before.
#'
#' @export
#' @seealso \link{qsf_filter} \link{qsm_message}
qsf_log = function(qsf)
{
  # 1. Extract messages from attributes
  messages <- lapply(qsf, function(x) attr(x, "message"))

  # 2. Identify non-empty messages
  keep_idx <- which(lengths(messages) > 0)

  if (length(keep_idx) == 0) return(list())

  ids <- names(qsf)

  # 3. Parse every message of every tree that has one, keeping track of the
  # tree index and treeID each message belongs to.
  rows <- list()
  for (i in keep_idx)
  {
    for (msg in messages[[i]])
    {
      parsed <- .qsf_parse_message(msg)
      rows[[length(rows) + 1]] <- list(
        index = i,
        treeID = suppressWarnings(as.integer(ids[i])),
        code = parsed$code,
        type = parsed$type,
        label = parsed$label,
        message = parsed$message
      )
    }
  }

  # 4. Group by code when available, otherwise by label (legacy behavior),
  # otherwise by the raw message as a last resort.
  keys <- vapply(rows, function(r)
  {
    if (!is.na(r$code)) return(r$code)
    if (!is.na(r$label)) return(r$label)
    return(r$message)
  }, character(1))

  unique_keys <- unique(keys)

  final_output <- lapply(unique_keys, function(key)
  {
    match_mask <- keys == key
    matches <- rows[match_mask]

    list(
      code    = matches[[1]]$code,
      type    = matches[[1]]$type,
      label   = matches[[1]]$label,
      message = matches[[1]]$message,
      index   = unique(vapply(matches, `[[`, integer(1), "index")),
      treeID  = unique(vapply(matches, `[[`, integer(1), "treeID"))
    )
  })

  names(final_output) <- unique_keys
  class(final_output) <- c("qsf_log", class(final_output))
  final_output
}

#' @rdname qsf_log
#' @param x A `qsf_log` object, as returned by \link{qsf_log}.
#' @param ... Unused (for S3 compatibility).
#' @export
print.qsf_log = function(x, ...)
{
  if (length(x) == 0)
  {
    cat("<qsf_log: no messages>\n")
    return(invisible(x))
  }

  for (nm in names(x))
  {
    entry <- x[[nm]]
    code  <- if (!is.na(entry$code)) entry$code else "--"
    type  <- if (!is.na(entry$type)) entry$type else "unknown"
    label <- if (!is.na(entry$label)) entry$label else nm
    cat(sprintf("[%s] %-8s %-30s n = %d trees\n", code, type, label, length(entry$treeID)))
  }

  invisible(x)
}

#' Filter a qsf object
#'
#' Subsets a `qsf` object, retaining only the trees that satisfy the
#' requested criteria. Filters compose with AND semantics: a tree is retained
#' only if it satisfies every criterion supplied.
#'
#' @param qsf A `qsf` object.
#' @param code character. One or more \link{qsf_log} codes (e.g. `c("W0", "E1")`).
#' By default, trees carrying any of these codes are **excluded**. Use `invert = TRUE`
#' to instead **keep only** the trees carrying one of these codes.
#' @param type character. One or more log types: `"warning"`, `"message"`, `"error"`.
#' Same semantics as `code`, and combined with `code` (if both given, an entry
#' must match both to be considered flagged).
#' @param invert logical. If `FALSE` (default) trees flagged by `code`/`type` are
#' dropped. If `TRUE` only the flagged trees are kept (useful for QA workflows).
#' @param dbh numeric(2). Keep only trees with a DBH (see \link{qsm_dbh}) inside
#' `c(min, max)`, in meters.
#' @param max_height numeric(2). Keep only trees with a height (see \link{qsm_height})
#' inside `c(min, max)`, in meters.
#' @param ... Unused. Reserved for additional filters in the future.
#'
#' @return A `qsf` object containing only the selected QSMs.
#'
#' @examples
#' \dontrun{
#' # Drop trees flagged with a warning code W0 or W3, and keep only merchantable trees
#' qsf_filter(qsf, code = c("W0", "W3"), dbh = c(0.09, Inf))
#'
#' # Keep only the trees flagged as errors, for inspection
#' qsf_filter(qsf, type = "error", invert = TRUE)
#' }
#'
#' @export
#' @seealso \link{qsf_log} \link{qsf_merchantable}
qsf_filter = function(qsf, ...)
{
  UseMethod("qsf_filter")
}

#' @export
#' @rdname qsf_filter
qsf_filter.qsf = function(qsf, code = NULL, type = NULL, invert = FALSE, dbh = NULL, max_height = NULL, ...)
{
  n <- length(qsf)
  keep <- rep(TRUE, n)

  # --- filter by log code and/or type ---
  if (!is.null(code) || !is.null(type))
  {
    log <- qsf_log(qsf)
    flagged <- rep(FALSE, n)

    for (entry in log)
    {
      match_code <- is.null(code) || (!is.na(entry$code) && entry$code %in% code)
      match_type <- is.null(type) || (!is.na(entry$type) && entry$type %in% type)

      if (match_code && match_type) flagged[entry$index] <- TRUE
    }

    keep <- keep & (if (invert) flagged else !flagged)
  }

  # --- filter by DBH range ---
  if (!is.null(dbh))
  {
    if (length(dbh) != 2) stop("`dbh` must be a numeric vector of length 2: c(min, max).", call. = FALSE)

    dbh_val <- vapply(qsf, function(x)
    {
      d <- tryCatch(qsm_dbh(x)$dbh, error = function(e) NA_real_)
      if (length(d) == 0) NA_real_ else d[1]
    }, numeric(1))

    keep <- keep & !is.na(dbh_val) & dbh_val >= dbh[1] & dbh_val <= dbh[2]
  }

  # --- filter by height range ---
  if (!is.null(max_height))
  {
    if (length(max_height) != 2) stop("`max_height` must be a numeric vector of length 2: c(min, max).", call. = FALSE)

    height_val <- qsm_height(qsf)

    keep <- keep & !is.na(height_val) & height_val >= max_height[1] & height_val <= max_height[2]
  }

  qsf[keep] |> as_qsf()
}

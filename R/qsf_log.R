#' QSF log
#'
#' Use qsf_log after \link{qsf} to get a structured summary of the log messages
#' recorded during the QSM reconstruction of every tree in a `qsf` object.
#'
#' @param qsf qsf
#'
#' @return
#' A named list, one element per distinct log entry. Elements are named by
#' their `code` (e.g. `"W0"`). Each element is a list with:
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
#' Messages recorded on a QSM start with a code `[<L><n>]`, where `<L>` is
#' one of `W` (warning), `M` (message) or `E` (error) and `<n>` is a numeric
#' identifier, immediately followed by `[Short label]` and the free text explanation,
#' e.g. `"[W0] [No wood point] This tree has no point labelled as wood"`. This is
#' fully backward compatible: messages recorded before this convention was introduced
#' have no code (`code = NA`) and are grouped by their short label instead, exactly
#' as before.
#' @examples
#' f <- system.file("extdata", "oak-plantation.qsf", package="arbor")
#' qsf <- qsf_read(f)
#' logs <- qsf_log(qsf)
#' print(logs)
#' logs$W2
#' logs$W2$treeID
#' logs$W3
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

#' Registry of known QSF log codes
#'
#' Internal registry mapping a log `code` (e.g. `"W0"`) to its default short
#' label. The leading letter of a code identifies its `type`: `W` = warning,
#' `M` = message, `E` = error. This registry is used as a fallback to label
#' coded messages that do not repeat their label in the message text, and to
#' document the codes produced by the underlying QSM reconstruction.
#'
#' @keywords internal
#' @noRd
.arbor_qsf_codes <- c(
  W0 = "No wood point",
  E0 = "QSM Generation Failed",
  W2 = "No valid measure",
  W3 = "Small tree allometry",
  M0 = "Broken tree",
  W5 = "Diameter anomaly"
)

#' Get the type (warning/message/error) associated with a log code
#' @param code character. A code such as `"W0"`.
#' @keywords internal
#' @noRd
.qsf_code_type <- function(code)
{
  if (is.na(code)) return(NA_character_)

  types <- c(
    W = "warning",
    M = "message",
    E = "error"
  )

  unname(types[substr(code, 1, 1)])
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
#' @noRd
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

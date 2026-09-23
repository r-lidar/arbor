#' Filter a qsf object
#'
#' `qsf_filter()` subsets a `qsf` object, retaining only the trees that
#' satisfy the requested criteria. Filters compose with AND semantics: a tree
#' is retained only if it satisfies every criterion supplied.\cr\cr
#' A handful of `qsf_filter_*()` wrappers cover the most common triage needs
#' on top of `qsf_filter()`: pulling out clean ("ok") trees, trees carrying
#' any warning, error or message, or (with `qsf_filter_flagged()`) anything
#' with a log entry at all. Each one still takes `...`, so it composes with
#' the other `qsf_filter()` criteria (`dbh`, `height`, ...) in a single
#' call. See details section for more in depth details
#' 
#' - `qsf_filter_ok()` retains QSMs with no log entry at all. They are assumed to be 
#'   good QSMs.
#' - `qsf_filter_sapling()` retains QSMs for which no measurements were performed
#'   according to the height allometry rules (see [Arbor Book](https://r-lidar.github.io/arbor_book/qsm.html#sec-tinytree)).
#'   Under Arbor's definition, they are considered saplings.
#' - `qsf_filter_nomeasure()` retains QSMs for which no measurements were found because
#'   the QSM algorithm failed to measure anything. However, these trees are not saplings.
#'   Saplings are different. For sapling the algortihm does not even try, on purpose.
#' - `qsf_filter_broken()` retains QSMs for which Arbor automatically detected
#'   that they correspond to broken trees.
#' - `qsf_filter_flagged()` retains QSMs with at least one log entry, of any kind.
#' - `qsf_filter_warnings()` retains QSMs carrying at least one warning.
#' - `qsf_filter_errors()` retains QSMs carrying at least one error.
#' - `qsf_filter_messages()` retains QSMs carrying at least one plain message.
#'
#' @param qsf A `qsf` object.
#' @param code character. One or more \link{qsf_log} codes (e.g. `c("W0", "E1")`).
#' By default, trees carrying any of these codes are **kept**. Use `invert = TRUE`
#' to instead **exclude** trees carrying these codes. In `qsf_filter_ok()`,
#' `qsf_filter_warnings()`, `qsf_filter_errors()` and `qsf_filter_messages()`,
#' `code` narrows the match to specific codes within that function's type;
#' left `NULL` (default) it matches any code of that type.
#' @param type character. One or more log types: `"warning"`, `"message"`, `"error"`.
#' Same semantics as `code`, and combined with `code` (if both given, an entry
#' must match both to be considered flagged). Only exposed on `qsf_filter()`
#' and `qsf_filter_flagged()`; the other wrappers fix `type` for you.
#' @param invert logical. If `FALSE` (default), trees matching `code`/`type` are
#' kept. If `TRUE`, matching trees are dropped. Only exposed on `qsf_filter()`.
#' @param flag logical. Filters on the mere *presence* of a log entry, regardless
#' of its code or type (independent of `code`/`type`, combined with them by AND
#' like every other filter here). `FALSE` keeps only trees with **no log entry
#' at all** (clean/"good" QSMs). `TRUE` keeps only trees with **at least one**
#' log entry, of any code or type. Default `NULL`: no filtering on this criterion.
#' This is what you want when you need "everything with nothing flagged" without
#' having to enumerate every possible `type`. Only exposed on `qsf_filter()`.
#' @param dbh numeric(2). Keep only trees with a DBH (see \link{qsm_dbh}) inside
#' `c(min, max)`, in meters.
#' @param height numeric(2). Keep only trees with a height (see \link{qsm_height})
#' inside `c(min, max)`, in meters.
#' @param ... In `qsf_filter()`: unused, reserved for additional filters in the
#' future. In `qsf_filter_*()` wrappers: passed on to `qsf_filter()` (e.g.
#' `dbh = c(0.2, Inf)`, `height = c(0, 30)`), so a wrapper can be chained
#' with the other filters in a single call.
#'
#' @return A `qsf` object containing only the selected QSMs.
#'
#' @examples
#' f <- system.file("extdata", "oak-plantation.qsf", package="arbor")
#' qsf <- qsf_read(f)
#'
#' # qsf_filter_*(): quick selectors
#' qsf_filter_ok(qsf)               # trees with no log entry at all
#' qsf_filter_sapling(qsf)          # saplings according to arbor's definition
#' qsf_filter_flagged(qsf)          # trees with at least one log entry, any kind
#' qsf_filter_warnings(qsf)         # trees carrying at least one warning
#' qsf_filter_warnings(qsf, "W2")   # ... specifically warning code W2
#' qsf_filter_messages(qsf)         # trees carrying at least one plain message
#'
#' # qsf_filter(): the general-purpose filter
#' # Keep only trees flagged with code W0 or W3
#' qsf_filter(qsf, code = c("W0", "W3"))
#'
#' # Exclude trees flagged with code W2 or W3, while keeping trees with DBH > 20 cm
#' qsf_filter(qsf, code = c("W2", "W3"), invert = TRUE, dbh = c(0.2, Inf))
#'
#' # Trees with no log entry at all ("good" QSMs)
#' qsf_filter(qsf, flag = FALSE)
#'
#' # Trees with at least one log entry, of any kind
#' qsf_filter(qsf, flag = TRUE)
#'
#' @export
#' @rdname qsf_filter
#' @seealso \link{qsf_log}
#' @md
qsf_filter <- function(qsf, code = NULL, type = NULL, invert = FALSE, flag = NULL, dbh = NULL, height = NULL, ...)
{
  n <- length(qsf)
  keep <- rep(TRUE, n)
  log <- NULL

  # --- filter by log code and/or type ---
  if (!is.null(code) || !is.null(type))
  {
    log <- qsf_log(qsf)
    flagged <- .qsf_log_match(log, code, type, n)

    # invert = FALSE keeps flagged trees; invert = TRUE drops flagged trees
    keep <- keep & (if (invert) !flagged else flagged)
  }

  # --- filter by mere presence/absence of ANY log entry, regardless of code/type ---
  if (!is.null(flag))
  {
    if (!is.logical(flag) || length(flag) != 1 || is.na(flag))
      stop("`flag` must be a single logical value (TRUE or FALSE).", call. = FALSE)

    if (is.null(log)) log <- qsf_log(qsf)
    any_log <- .qsf_log_match(log, NULL, NULL, n)

    keep <- keep & (if (flag) any_log else !any_log)
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
  if (!is.null(height))
  {
    if (length(height) != 2) stop("`height` must be a numeric vector of length 2: c(min, max).", call. = FALSE)

    height_val <- qsm_height(qsf)

    keep <- keep & !is.na(height_val) & height_val >= height[1] & height_val <= height[2]
  }

  qsf[keep]
}

#' @export
#' @rdname qsf_filter
qsf_filter_ok <- function(qsf, ...)
{
  qsf_filter(qsf, flag = FALSE, ...)
}

#' @export
#' @rdname qsf_filter
qsf_filter_sapling <- function(qsf, ...)
{
  qsf_filter(qsf, code = "W3", ...)
}

#' @export
#' @rdname qsf_filter
qsf_filter_nomeasure <- function(qsf, ...)
{
  qsf_filter(qsf, code = "W2", ...)
}

#' @export
#' @rdname qsf_filter
qsf_filter_broken <- function(qsf, ...)
{
  qsf_filter(qsf, code = "W4")
}

#' @export
#' @rdname qsf_filter
qsf_filter_flagged <- function(qsf, code = NULL, type = NULL, ...)
{
  if (is.null(code) && is.null(type))
    qsf_filter(qsf, flag = TRUE, ...)
  else
    qsf_filter(qsf, code = code, type = type, ...)
}

#' @export
#' @rdname qsf_filter
qsf_filter_warnings <- function(qsf, code = NULL, ...)
{
  qsf_filter(qsf, type = "warning", code = code, ...)
}

#' @export
#' @rdname qsf_filter
qsf_filter_errors <- function(qsf, code = NULL, ...)
{
  qsf_filter(qsf, type = "error", code = code, ...)
}

#' @export
#' @rdname qsf_filter
qsf_filter_messages <- function(qsf, code = NULL, ...)
{
  qsf_filter(qsf, type = "message", code = code, ...)
}

# Internal helper: logical vector, length n, TRUE where a log entry matches
# both `code` and `type` (NULL = match anything). Shared by the code/type
# block and the `flag` block of qsf_filter() so the log is only walked once
# per criterion and the matching logic lives in exactly one place.
.qsf_log_match <- function(log, code, type, n)
{
  matched <- rep(FALSE, n)

  for (entry in log)
  {
    match_code <- is.null(code) || (!is.na(entry$code) && entry$code %in% code)
    match_type <- is.null(type) || (!is.na(entry$type) && entry$type %in% type)

    if (match_code && match_type && !is.null(entry$index))
    {
      idx <- entry$index[!is.na(entry$index) & entry$index >= 1 & entry$index <= n]
      matched[idx] <- TRUE
    }
  }

  matched
}
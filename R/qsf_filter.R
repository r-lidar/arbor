#' Filter a qsf object
#'
#' Subsets a `qsf` object, retaining only the trees that satisfy the
#' requested criteria. Filters compose with AND semantics: a tree is retained
#' only if it satisfies every criterion supplied.
#'
#' @param qsf A `qsf` object.
#' @param code character. One or more \link{qsf_log} codes (e.g. `c("W0", "E1")`).
#' By default, trees carrying any of these codes are **kept**. Use `invert = TRUE`
#' to instead **exclude** trees carrying these codes.
#' @param type character. One or more log types: `"warning"`, `"message"`, `"error"`.
#' Same semantics as `code`, and combined with `code` (if both given, an entry
#' must match both to be considered flagged).
#' @param invert logical. If `FALSE` (default), trees matching `code`/`type` are
#' kept. If `TRUE`, matching trees are dropped.
#' @param dbh numeric(2). Keep only trees with a DBH (see \link{qsm_dbh}) inside
#' `c(min, max)`, in meters.
#' @param max_height numeric(2). Keep only trees with a height (see \link{qsm_height})
#' inside `c(min, max)`, in meters.
#' @param ... Unused. Reserved for additional filters in the future.
#'
#' @return A `qsf` object containing only the selected QSMs.
#'
#' @examples
#' f <- system.file("extdata", "oak-plantation.qsf", package="arbor")
#' qsf <- qsf_read(f)
#' 
#' # Keep only trees flagged with code W2
#' qsf_select(qsf, code = c("W0", "W3"))
#'
#' # Exclude trees flagged with code W2 or W3, while keeping trees with DBH > 20 cm
#' qsf_select(qsf, code = c("W2", "W3"), invert = TRUE, dbh = c(0.2, Inf))
#' @export
#' @seealso \link{qsf_log} \link{qsf_merchantable}
qsf_select <- function(qsf, code = NULL, type = NULL, invert = FALSE, dbh = NULL, max_height = NULL, ...)
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

      if (match_code && match_type && !is.null(entry$index))
      {
        idx <- entry$index[!is.na(entry$index) & entry$index >= 1 & entry$index <= n]
        flagged[idx] <- TRUE
      }
    }

    # invert = FALSE keeps flagged trees; invert = TRUE drops flagged trees
    keep <- keep & (if (invert) !flagged else flagged)
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
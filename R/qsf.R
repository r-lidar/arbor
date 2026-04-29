#' Quantitative Structural Forest
#'
#' Batch processing of QSM models with parallel execution.
#'
#' @param las A point cloud with semantic and instance segmentation computed
#' @param params list See \link{parameters}.
#' @return A qsf object
#'
#' @export
#' @seealso  \link{qsm}
qsf <- function(las, params = default_arbor_parameters)
{
  #params <- evaluate_penalty(params)
  res <- qsf_cpp(las@data, params)
  for(i in seq_along(res)) res[[i]] <- qsm_finalize(res[[i]])
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

#' QSF log
#'
#' Use qsf_log() after qsf() to get the logs
#'
#' @param qsf qsf
#' @export
qsf_log = function(qsf)
{
  messages = lapply(qsf, function(x) attr(x, "message"))

  # Identify which elements are NOT empty
  keep_idx <- which(sapply(messages, length) > 0)

  # Subset messages and store their original positions
  clean_list <- messages[keep_idx]

  # Extract the tags to use as grouping keys
  warn_tags <- sapply(clean_list, function(x) {
    regmatches(x, regexpr("\\[WARN \\d+\\]", x))
  })

  # Create the structured list
  unique_tags <- unique(unlist(warn_tags))

  final_output <- lapply(unique_tags, function(tag)
  {
    match_mask <- warn_tags == tag

    matches <- clean_list[match_mask]
    original_indices <- keep_idx[match_mask]

    # Get the first message and "templatize" it
    # This replaces digits/decimals with 'x' to make it generic
    raw_msg <- matches[[1]]
    generic_msg <- gsub("\\d+\\.\\d+|\\d+", "x", raw_msg)
    generic_msg <- trimws(gsub("\\[WARN x]", "", generic_msg))

    list(
      message = generic_msg,
      index = unname(original_indices),
      treeID = as.integer(names(matches))
    )
  })

  # Name the list elements by their tag
  names(final_output) <- unique_tags
  final_output
}

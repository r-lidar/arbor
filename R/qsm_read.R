#' Read QSM Data from File
#'
#' Loads QSM segment data from a CSV file and unify naming convention. So the function reads Computree
#' TreeQSM, arbor in the same format.
#'
#' @param x A character string specifying the file path.
#' @return A data frame of QSM segment data.
#' @export
#' @md
qsm_read = function(x)
{
  qsm = data.table::fread(x)
  qsm = unify_names(qsm)
  name = tools::file_path_sans_ext(basename(x))
  attr(qsm, "ID") = name
  return(qsm)
}

unify_names <- function(qsm)
{
  # Define the renaming mapping: old = new
  name_map <- c(
    "branchOrder" = "branch_order",
    "parent_segment_ID" = "parent_ID",
    "segment_ID" = "cyl_ID"
  )

  # Replace column names only if they exist in the data.frame
  current_names <- names(qsm)
  new_names <- ifelse(current_names %in% names(name_map), name_map[current_names], current_names)
  names(qsm) <- new_names
  qsm
}



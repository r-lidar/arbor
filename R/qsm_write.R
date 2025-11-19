#' Write a QSM to File
#'
#' Writes a Quantitative Structure Model (QSM) to a file in either OBJ, PLY or CSV format,
#' based on the file extension.
#'
#' Supported formats:
#' * `.ply` or `.obj`: writes the QSM as a mesh file
#' * `.csv` or `.txt`: writes the QSM as a table
#'
#' @param qsm A QSM object to be written.
#' @param file A string giving the path to the output file. The file extension determines
#'   the format (e.g., ".ply", "obj", ".csv", ".txt").
#'
#' @export
#'
#' @examples
#' \dontrun{
#' qsm_write(qsm, "tree.ply")
#' qsm_write(qsm, "tree.csv")
#' }
#' @export
#' @md
qsm_write = function(qsm, file)
{
  file = normalizePath(file, mustWork = FALSE)
  ext = tools::file_ext(file)
  if (ext %in% c("ply", "obj"))
    qsm_write_cpp(qsm, file)
  else if (ext %in% c("csv", "txt"))
    qsm_write_table(qsm, file)
  else
    stop("format not supported")

  return(invisible(TRUE))
}

qsm_write_table = function(qsm, file)
{
  data.table::fwrite(qsm, file)
}


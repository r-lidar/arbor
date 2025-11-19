#' Write a QSM to File
#'
#' Writes a Quantitative Structure Model (QSM) to a file in either OBJ, PLY or CSV format,
#' based on the file extension.
#'
#' Supported formats:
#' * `.ply` or `.obj` or `.stl`: writes the QSM as a mesh file
#' * `.csv` or `.txt`: writes the QSM as a table
#'
#' @param qsm A QSM object to be written.
#' @param file A string giving the path to the output file. The file extension determines
#'   the format (e.g., ".ply", "obj", ".csv", ".txt", ".stl").
#' @param binary Boolean. Used if the format supports ASCII or binary
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
qsm_write = function(qsm, file, binary = TRUE)
{
  file = normalizePath(file, mustWork = FALSE)
  qsm_write_cpp(qsm, file, binary)
  return(invisible(TRUE))
}

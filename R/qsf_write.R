#' Write a QSF to Files
#'
#' Writes a Quantitative Structure Forest (QSF) to files in either OBJ, PLY or CSV format,
#' based on the file extension. Each QSM is written in its own file.
#'
#' Supported formats:
#' * `.ply` or `.obj` or `.stl`: writes the QSM as a mesh file
#' * `.csv` or `.txt`: writes the QSM as a table
#'
#' @param qsf A QSF object to be written.
#' @param dir A string giving the director to the output files.
#' @param format the format (e.g., ".ply", "obj", ".csv", ".txt", ".stl").
#' @param binary Boolean. Used if the format supports ASCII or binary
#'
#' @export
#'
#' @examples
#' \dontrun{
#' qsf_write(qsf, "directory/csv", format = "csv")
#' qsf_write(qsf, "directory/obj", format = "obj")
#' }
#' @export
#' @md
qsf_write = function(qsf, dir, format = "csv", binary = TRUE)
{
  file = normalizePath(dir, mustWork = FALSE)
  qsf_write_cpp(qsf, file, binary)
  return(invisible(TRUE))
}

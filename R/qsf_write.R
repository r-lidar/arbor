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
#' @param formats the format (e.g., ".ply", "obj", ".csv", ".txt", ".stl").
#' @param binary Boolean. Used if the format supports ASCII or binary
#'
#' @export
#' @export
#' @md
qsf_write = function(qsf, dir, formats = c("csv", "obj"), binary = TRUE)
{
  dir = normalizePath(dir, mustWork = FALSE)
  for (format in formats) qsf_write_cpp(qsf, dir, format, binary)
  return(invisible(TRUE))
}

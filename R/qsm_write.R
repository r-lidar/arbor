#' Write a QSM to File
#'
#' Writes a Quantitative Structure Model (QSM) to a file in either PLY or CSV format,
#' based on the file extension.
#'
#' Supported formats:
#' - `.ply`: writes the QSM as a PLY mesh file
#' - `.obj`: writes the QSM as a OBJ mesh file
#' - `.csv` or `.txt`: writes the QSM as a table
#'
#' @param qsm A QSM object to be written.
#' @param file A string giving the path to the output file. The file extension determines
#'   the format (e.g., ".ply", ".csv", ".txt").
#'
#' @export
#'
#' @examples
#' \dontrun{
#' qsm_write(qsm, "tree.ply")
#' qsm_write(qsm, "tree.csv")
#' }
#' @export
qsm_write = function(qsm, file)
{
  ext = tools::file_ext(file)
  if (ext %in% c("ply"))
    write_qsm_ply(qsm, file)
  else if (ext %in% c("csv", "txt"))
    write_qsm_table(qsm, file)
  else if (ext %in% c("obj"))
    write_qsm_obj(qsm, file)
  else
    stop("format not supported")

  return(invisible(TRUE))
}


write_qsm_table = function(qsm, file)
{
  data.table::fwrite(qsm, file)
}


write_qsm_ply <- function(df, file, resolution = 16L)
{
  normalize <- function(v) v / sqrt(sum(v^2))
  cross <- function(a, b) c(
    a[2]*b[3] - a[3]*b[2],
    a[3]*b[1] - a[1]*b[3],
    a[1]*b[2] - a[2]*b[1]
  )

  vertices <- list()
  faces <- list()
  v_offset <- 0L

  for (i in seq_len(nrow(df)))
  {
    row <- df[i, ]
    p0 <- as.numeric(c(row$startX, row$startY, row$startZ))
    p1 <- as.numeric(c(row$endX, row$endY, row$endZ))
    radius <- row$radius

    axis <- p1 - p0
    if (all(axis == 0)) next  # Skip degenerate cylinders

    height <- sqrt(sum(axis^2))
    dir <- normalize(axis)
    z_axis <- c(0, 0, 1)

    # Build orthonormal basis
    ortho1 <- if (all(abs(dir - z_axis) < 1e-6)) c(1, 0, 0) else normalize(cross(z_axis, dir))
    ortho2 <- normalize(cross(dir, ortho1))

    # Create base circle
    angles <- seq(0, 2 * pi, length.out = resolution + 1L)[-1L]
    circle <- t(sapply(angles, function(theta)
    {
      cos(theta) * ortho1 + sin(theta) * ortho2
    })) * radius

    # Compute bottom and top vertices
    bottom <- sweep(circle, 2, p0, FUN = "+")
    top <- sweep(circle, 2, p1, FUN = "+")
    verts <- rbind(bottom, top)
    vertices[[i]] <- verts

    # Create faces
    for (j in seq_len(resolution))
    {
      a <- v_offset + j
      b <- v_offset + if (j < resolution) j + 1L else 1L
      c <- a + resolution
      d <- b + resolution

      # Each quad is two triangles
      faces[[length(faces) + 1L]] <- c(a - 1L, b - 1L, d - 1L)
      faces[[length(faces) + 1L]] <- c(a - 1L, d - 1L, c - 1L)
    }

    v_offset <- v_offset + 2L * resolution
  }

  vertices <- do.call(rbind, vertices)
  faces <- do.call(rbind, faces)

  vertices <- round(vertices, 3)

  # Write to PLY (ASCII)
  con <- file(file, "w")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "ply",
    "format ascii 1.0",
    paste("element vertex", nrow(vertices)),
    "property double x",
    "property double y",
    "property double z",
    paste("element face", nrow(faces)),
    "property list uchar int vertex_indices",
    "end_header"
  ), con)

  # Use sprintf to format with exactly 3 decimals
  writeLines(apply(vertices, 1, function(x) sprintf("%.3f %.3f %.3f", x[1], x[2], x[3])), con)
  utils::write.table(cbind(3L, faces), file = con, row.names = FALSE, col.names = FALSE)
}

write_qsm_obj <- function(df, file, resolution = 16L)
{
  normalize <- function(v) v / sqrt(sum(v^2))
  cross <- function(a, b) c(
    a[2]*b[3] - a[3]*b[2],
    a[3]*b[1] - a[1]*b[3],
    a[1]*b[2] - a[2]*b[1]
  )

  vertices <- list()
  faces <- list()
  v_offset <- 0L

  for (i in seq_len(nrow(df)))
  {
    row <- df[i, ]
    p0 <- as.numeric(c(row$startX, row$startY, row$startZ))
    p1 <- as.numeric(c(row$endX, row$endY, row$endZ))
    radius <- row$radius

    axis <- p1 - p0
    if (all(axis == 0)) next

    dir <- normalize(axis)
    z_axis <- c(0, 0, 1)

    ortho1 <- if (all(abs(dir - z_axis) < 1e-6)) c(1, 0, 0) else normalize(cross(z_axis, dir))
    ortho2 <- normalize(cross(dir, ortho1))

    angles <- seq(0, 2 * pi, length.out = resolution + 1L)[-1L]
    circle <- t(sapply(angles, function(theta)
    {
      cos(theta) * ortho1 + sin(theta) * ortho2
    })) * radius

    bottom <- sweep(circle, 2, p0, FUN = "+")
    top <- sweep(circle, 2, p1, FUN = "+")
    verts <- rbind(bottom, top)
    vertices[[length(vertices) + 1L]] <- verts

    for (j in seq_len(resolution))
    {
      a <- v_offset + j
      b <- v_offset + if (j < resolution) j + 1L else 1L
      c <- a + resolution
      d <- b + resolution

      # OBJ is 1-based
      faces[[length(faces) + 1L]] <- c(a, b, d)
      faces[[length(faces) + 1L]] <- c(a, d, c)
    }

    v_offset <- v_offset + 2L * resolution
  }

  # Combine all vertex and face data
  vertices <- do.call(rbind, vertices)
  faces <- do.call(rbind, faces)

  con <- file(file, "w")
  on.exit(close(con), add = TRUE)

  # Write vertices
  for (i in seq_len(nrow(vertices))) {
    v <- vertices[i, ]
    writeLines(sprintf("v %.3f %.3f %.3f", v[1], v[2], v[3]), con)
  }

  # Write faces
  for (i in seq_len(nrow(faces))) {
    f <- faces[i, ]
    writeLines(sprintf("f %d %d %d", f[1], f[2], f[3]), con)
  }
}



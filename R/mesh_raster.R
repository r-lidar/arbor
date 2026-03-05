write_raster_to_obj <- function(r, filename, z_scale = 1)
{
  nr <- terra::nrow(r)
  nc <- terra::ncol(r)

  # --- STEP 1: Filter Vertices ---
  # Get indices of all non-NA cells
  vals <- terra::values(r, mat = FALSE)
  valid_idx <- which(!is.na(vals))

  if (length(valid_idx) == 0) {
    stop("No valid (non-NA) cells in raster.")
  }

  # Extract only valid coordinates and z-values
  coords <- terra::xyFromCell(r, valid_idx)
  z_vals <- vals[valid_idx]

  vertices <- cbind(coords, z_vals * z_scale)

  # --- STEP 2: Create Index Map ---
  # This map translates: Original_Cell_Index -> New_Vertex_ID
  id_map <- rep(NA_integer_, terra::ncell(r))
  id_map[valid_idx] <- seq_along(valid_idx)

  # --- STEP 3: Build Faces ---
  # Iterate through grid cells to find valid quads
  face_list <- vector("list", (nc - 1) * (nr - 1))
  face_count <- 0

  for (j in 1:(nc - 1)) {
    for (i in 1:(nr - 1)) {
      # Cell indices for the 4 corners of a quad
      # Terra uses row-major order: cell_index = (row - 1) * ncol + col

      # Bottom-left, top-left, top-right, bottom-right
      c_bl <- (i - 1) * nc + j       # bottom-left (row i, col j)
      c_tl <- i * nc + j             # top-left (row i+1, col j)
      c_tr <- i * nc + (j + 1)       # top-right (row i+1, col j+1)
      c_br <- (i - 1) * nc + (j + 1) # bottom-right (row i, col j+1)

      # Get remapped vertex IDs
      v_ids <- id_map[c(c_bl, c_tl, c_tr, c_br)]

      # Only create face if ALL 4 corners are valid (non-NA)
      if (!anyNA(v_ids)) {
        face_count <- face_count + 1
        face_list[[face_count]] <- v_ids
      }
    }
  }

  # --- STEP 4: Write File ---
  if (face_count == 0) {
    stop("No valid faces to write. All quads contain at least one NA cell.")
  }

  # Trim unused list elements
  face_list <- face_list[1:face_count]
  faces_final <- do.call(rbind, face_list)

  # Write OBJ file
  con <- file(filename, "w")
  on.exit(close(con), add = TRUE)

  cat("# OBJ Export: NA cells excluded\n", file = con)
  cat(sprintf("# Generated: %s\n", Sys.time()), file = con)
  cat(sprintf("# Valid vertices: %d | Valid faces: %d\n", nrow(vertices), nrow(faces_final)), file = con)

  # Write vertices
  utils::write.table(cbind("v", vertices), file = con, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

  # Write faces (OBJ uses 1-based indexing, which we already have)
  utils::write.table(cbind("f", faces_final), file = con, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

  message(sprintf("File saved: %s", filename))
  message(sprintf("Vertices: %d | Faces: %d | Coverage: %.1f%%", nrow(vertices), nrow(faces_final), 100 * length(valid_idx) / terra::ncell(r)))

  invisible(list(vertices = nrow(vertices), faces = nrow(faces_final)))
}

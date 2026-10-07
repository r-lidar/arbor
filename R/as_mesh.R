as_mesh <- function(x, ...) {
  UseMethod("as_mesh")
}

as_mesh.qsm <- function(x, offset = c(0, 0), color = "cyl_ID", pal = "auto", precomputed_colors = NULL, ...) {
  if (nrow(x) == 0) return(NULL)
  
  # Apply XY translation to coordinates before generating the mesh
  if (!is.null(offset)) {
    x$startX <- x$startX - offset[1]
    x$startY <- x$startY - offset[2]
    x$endX   <- x$endX - offset[1]
    x$endY   <- x$endY - offset[2]
  }

  pal <- get_qsm_pal(color, pal)

  # Color mapping
  if (!is.null(precomputed_colors)) {
    colors_mapped <- precomputed_colors
  } else {
    if (color %in% names(x)) {
      is_cat <- color %in% c("branch_order", "quality")
      colors_mapped <- map_qsm_colors(x[[color]], pal, is_cat)
    } else {
      colors_mapped <- rep("black", nrow(x))
    }
  }

  mesh_data <- qsm_mesh_cpp(x, 16)

  # Align colors with mesh vertices
  id <- match(mesh_data$NodeID, x$cyl_ID)
  Final_Colors <- colors_mapped[id]

  mesh <- rgl::qmesh3d(
    vertices = mesh_data$vertices,
    indices  = mesh_data$indices,
    homogeneous = FALSE
  )

  mesh$material$color <- Final_Colors
  mesh$material$specular <- "black"
  mesh$material$shininess <- 0

  return(mesh)
}

as_mesh.qsf <- function(x, offset = NULL, color = "branch_order", pal = "auto", ...) {
  # If no offset is provided, calculate the global minimums across all QSMs
  if (is.null(offset)) {
    tx <- min(sapply(x, function(qsm) min(qsm$startX)))
    ty <- min(sapply(x, function(qsm) min(qsm$startY)))
    offset <- c(tx, ty)
  }
  
  # Apply as_mesh to each QSM in the forest
  meshes <- lapply(x, as_mesh.qsm, offset = offset, color = color, pal = pal, ...)
  names(meshes) <- names(x)
  meshes <- Filter(Negate(is.null), meshes)
  
  # Attach the calculated/used offset as an attribute so plot.qsf can retrieve it
  attr(meshes, "offset") <- offset
  
  return(meshes)
}

get_qsm_pal <- function(color, pal) {
  if (identical(pal, "auto")) {
    default_pal <- c("blue", "green", "yellow", "orange", "red")
    branch_order_pal <- c("#552203","#ad3a01","#e59b16","#fdd63b","#9acd56", "#59ad4b","#1e490e")
    quality_pal <- c("blue", "green", "yellow", "orange", "red")
    
    pal <- switch(color, 
                  branch_order = branch_order_pal, 
                  quality = quality_pal, 
                  default_pal)
  }
  return(pal)
}
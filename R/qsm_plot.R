#' Plot a Cylinder in 3D Space
#'
#' Creates and plots a 3D cylinder between two points with a specified radius and color.
#'
#' @param start A numeric vector of length 3 representing the starting point of the cylinder.
#' @param end A numeric vector of length 3 representing the ending point of the cylinder.
#' @param radius A numeric value representing the radius of the cylinder.
#' @param color A character string for the cylinder color (default is "blue").
#' @export
#' @md
plot_cylinder <- function(start, end, radius, color="blue")
{
  # Create a cylinder object from start to end with specified radius
  cyl <- rgl::cylinder3d(
    radius = radius,
    sides = 10,
    # Define the two end points
    v1 = start,
    v2 = end
  )

  # Add the cylinder to the 3D plot
  rgl::shade3d(cyl, color = color, alpha = 0.5)
}

#' Plot QSM Data in 3D
#'
#' Plots the QSM data as a series of connected cylinders, with color-coding based on branch attributes.
#'
#' @param qsm A data frame containing QSM data with segment attributes.
#' @param add A numeric vector for translation offsets (default is NULL).
#' @param sides Number of sides for each cylinder (default is 16).
#' @param color The attribute for color mapping (default is "branchOrder").
#' @param name description
#' @param skeleton,cylinder boolean. plot the skeleton or the cylinders or both.
#' @export
#' @md
plot_qsm = function(qsm, add = NULL, sides = 16, color = "cyl_ID", skeleton = TRUE, cylinder = TRUE)
{
  color_palette <- grDevices::colorRampPalette(c("blue", "green", "yellow", "orange", "red"))
  colattr = qsm[[color]]

  if (!"radius" %in% names(qsm))
    cylinder = FALSE

  if (color %in% names(qsm))
  {
    if (color == "branch_order") {
      colattr = colattr
    } else if (is.logical(colattr)) {
      colattr = colattr+1
    } else {
      colattr = findInterval(colattr, seq(min(colattr), max(colattr), length.out = 20))
    }
    colors <- color_palette(max(colattr))
  }
  else
  {
    colattr = rep(1, nrow(qsm))
    colors = "black"
  }

  if (!"translateX" %in% names(qsm))
      qsm$translateX = min(qsm$startX)
  if (!"translateY" %in% names(qsm))
      qsm$translateY = min(qsm$startY)
  if (!"translateZ" %in% names(qsm))
    qsm$translateZ = 0


  if (!is.null(add))
  {
    tx = add[1]
    ty = add[2]
    tz = 0
  }
  else
  {
    x = -c(0,0)
    tx = qsm$translateX[1]
    ty = qsm$translateY[1]
    tz = 0
    rgl::open3d()
  }

  # Calculate axis limits
  z_limits <- as.integer(range(c(qsm$startZ-qsm$translateZ, qsm$endZ-qsm$translateZ)) + c(-1, 1))

  if (cylinder)
  {
    cyls = lapply(1:nrow(qsm), function(i)
    {
      start <- c(qsm$startX[i] - tx, qsm$startY[i] - ty, qsm$startZ[i] - tz)
      end <- c(qsm$endX[i] - tx, qsm$endY[i] - ty, qsm$endZ[i]- tz)
      radius <- qsm$radius[i]
      color <- colors[colattr[i]]

      m = rbind(start, end)
      cyl = rgl::cylinder3d(m, radius, closed = -1, sides = sides)
      cyl$material$color <- color
      cyl
    })
  }

  if (skeleton)
  {
    sktln = lapply(1:nrow(qsm), function(i)
    {
      start <- c(qsm$startX[i] - tx, qsm$startY[i] - ty, qsm$startZ[i] - tz)
      end <- c(qsm$endX[i] - tx, qsm$endY[i] - ty, qsm$endZ[i]- tz)

      m = rbind(start, end)
    })
    sktln = do.call(rbind, sktln)
  }

  if (skeleton)
  {
    qq = qsm[,1:3]
    qq$startX = qq$startX - tx
    qq$startY = qq$startY - ty
    qq$startZ = qq$startZ - tz
    rgl::segments3d(sktln, col = rep(colors[colattr], each = 2))
    rgl::points3d(qq, col = colors[colattr])
  }
  if (cylinder)
  {
    mesh = do.call(merge, cyls)
    rgl::shade3d(mesh)
  }

  # Define custom ticks for the z-axis
  z_ticks <- seq(z_limits[1], z_limits[2], by = 0.5)  # Finer ticks for z-axis
  z_labels <- as.character(z_ticks) # Labels as character strings

  # Manually add finer labels to the z-axis
  rgl::axis3d("z", at = z_ticks, labels = z_labels, col = "black")
  rgl::axis3d("x", col = "black")
  rgl::axis3d("y", col = "black")

  lidR:::.pan3d(2)

  return(invisible(c(tx, ty)))
}

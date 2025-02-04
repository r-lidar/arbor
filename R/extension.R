#' Extent the trees to the DTM
#'
#' Extent the trees to the DTM. It fits circles using a RANSAC approach to measure the bottom diameter
#' of each tree. In practice fitting one circle on a low slice of a stem works only in easy contexts with
#' vertical and perfectly segmented wood and no noise or mistake. Actaully the method fit numerous
#' circles at different height and with different thicknesses to be more robust and avoid edge cases
#' and pit-falls.
#'
#' @param trees LAS object from lidR
#' @param dtm SpatRaster from terra. The DTM.
#' @param max_diameter The maximum diameter of a tree. This only serve to protect against bad fitting
#' with small trees where RANSAC fitting could easy find too large diameters. Input a number that
#' correspond realistically to the size of a biggest tree in your context.
#' @param ... unused. Only useful to separate important parameters from those that should only be change
#' in very specific contexts
#' @param name description
#' @param max_height make slices for the bottom of the tree (i.e. 50 cm if ou sliced at 50 cm) up to
#' +max_height above that point (i.e. 1.5 m in this example).
#' @param min_slice_thickness,max_slice_thickness the algorithm fits circles on many slices with various
#' thicknesses. This control the thicknesses limits
#' @param extra_height the trees can be prolongated a below the ground
#' @param debug boolean. Some debugging display.
#' @export
tree_extensions = function(trees, dtm, max_diameter = 0.5, ...., max_height = 1, min_slice_thickness = 0.1, max_slice_thickness = 0.6, extra_height = 0, debug = FALSE)
{
  max_radius = max_diameter/2

  # Generate multiple slices
  ranges = expand.grid(bottom = seq(0.1, max_height, 0.05), top = seq(0.1, max_height, 0.05))
  ranges = ranges[ranges$top - ranges$bottom >= min_slice_thickness,]
  ranges = ranges[ranges$top - ranges$bottom <= max_slice_thickness,]

  cat("RANSAC fitting on", nrow(ranges), "slices per trees\n")

  extensions = list()
  for (id in unique(trees$treeID))
  {
    cat("Tree", id, ": ")
    tt = lidR::filter_poi(trees, treeID == id, foliage == FALSE)
    min_hag = min(tt$hag)
    tt = lidR::filter_poi(tt, hag < min_hag + max_height)

    #xyz = sf::st_coordinates(tt)
    #pca <- prcomp(xyz, center = TRUE, scale. = FALSE)
    #main_axis <- pca$rotation[, 1]  # First principal component

    # Compute the rotation matrix
    #rotation_matrix <- align_to_z(main_axis)

    # Apply the rotation to the point cloud
    #rotated_xyz = xyz %*% t(rotation_matrix)
    #rotated_xyz <- scale(rotated_xyz, center = TRUE, scale = FALSE)
    #ttt <- as.data.frame(rotated_xyz)
    #names(ttt) = c("X", "Y", "Z")
    #ttt$clusterID = tt$clusterID
    #ttt = LAS(ttt)


    # Plot the point cloud
    #rgl::plot3d(centered_xyz, col = "blue", size = 2)
    #rgl::points3d(rotated_xyz, col = "red", size = 2)
    ##rgl::arrow3d(p0 = c(0,0,0), p1 =  main_axis, type = "lines",  col = "red", length = 2)

    circles = apply(ranges, 1, function(x)
    {
      bottom = lidR::filter_poi(tt, hag  >= x[1], hag <= x[2])
      if (lidR::is.empty(bottom)) return(NULL)

      bottom$Z = bottom$Z * 0.01
      bottom = lidR::connected_components(bottom, 0.01, 5)
      ids = table(bottom$clusterID)
      ids = as.numeric(names(ids[which.max(ids)]))
      bottom = bottom[bottom$clusterID == ids]

      if (lidR::npoints(bottom) < 10) return(NULL)

      bottom$Z = bottom$Z * 100
      circle = lidR::fit_circle(bottom)
      inliner = length(circle$inliers)/npoints(bottom) * 100

      if (debug)
      {
        plot(sf::st_coordinates(bottom), asp = 1, main = id)
        graphics::symbols(circle$center_x, circle$center_y, circles = circle$radius,  add = TRUE, fg = "red", inches = FALSE)
        graphics::symbols(circle$center_x, circle$center_y, circles = circle$radius+0.01,  add = TRUE, fg = "red", lty=3, inches = FALSE)
        graphics::symbols(circle$center_x, circle$center_y, circles = circle$radius-0.01,  add = TRUE, fg = "red", lty= 3, inches = FALSE)
        graphics::mtext(paste0("Radius = ", round(circle$radius, 3), " inliner = ", round(inliner), "% sector ", circle$angle_range, " deg"))
      }
      circle$pinlier = inliner
      circle$inliers = NULL
      as.data.frame(circle)
    })

    if (is.null(circles)) next

    circles = do.call(rbind, circles)

    valid = circles$radius < max_radius & circles$angle_range > 90
    n = sum(valid)
    if (n == 0) valid = circles$radius < max_radius
    n = sum(valid)
    if (n == 0)
    {
      cat("no valid circle detected\n")
      next
    }

    circles = circles[valid,]
    circles = circles[rev(order(circles$pinlier)),]
    n = min(c(5, n))

    if (n == 0) next

    circles = circles[1:n,]
    x = stats::median(circles$center_x)
    y = stats::median(circles$center_y)
    z = min(circles$z)+0.05
    r = stats::median((circles$radius))
    circle = list(center_x = x, center_y = y, z = z, radius = r)

    loc = matrix(c(x,y), ncol = 2)
    zgnd = terra::extract(dtm, loc, method = "bilinear")
    hcyl = as.numeric(z-zgnd+extra_height)

    if (circle$radius > max_radius)
    {
      cat("radius =", round(circle$radius*100, 0), "cm is to big. No extension.")
    }
    else(circle$radius < max_radius)
    {
      extension = generate_cylinder_points(circle, height = hcyl)
      #extension = as.matrix(extension)
      #extension <- extension %*% rotation_matrix
      #extension = as.data.frame(extension)
      #names(extension) = c("X", 'Y', "Z")
      lidR::quantize(extension$X, tt@header[["X scale factor"]], tt@header[["X offset"]])
      lidR::quantize(extension$Y, tt@header[["Y scale factor"]], tt@header[["Y offset"]])
      lidR::quantize(extension$Z, tt@header[["Z scale factor"]], tt@header[["Z offset"]])
      extension$treeID = id
      extensions[[as.character(id)]] = extension

      cat("radius =", round(circle$radius*100, ), "cm")
    }

    cat("\n")
    #x = plot(tt)
    #plot(LAS(extension), add = x)
  }

  extensions = do.call(rbind, extensions)
  extensions$anisotropy = 1
  extensions$wood = TRUE
  data.table::setDT(extensions)

  extensions = lidR::LAS(extensions)

  return(extensions)
}

#' @rdname tree_extension
#' @param tree LAS object
#' @param extension object returned by \link{tree_extension}
#' @param fill_value default value to fill missing attributes when merging
#' @export
weld_extension <- function(trees, extensions, fill_value = 0L)
{
  ..all_cols <- NULL

  df1 = trees@data
  df2 = extensions@data

  all_cols <- union(names(df1), names(df2))

  # Add missing columns with the fill value
  for (col in setdiff(all_cols, names(df1))) {
    df1[[col]] <- fill_value
  }
  for (col in setdiff(all_cols, names(df2))) {
    df2[[col]] <- fill_value
  }

  # Ensure column order matches
  df1 <- df1[, ..all_cols]
  df2 <- df2[, ..all_cols]

  # Combine rows
  trees@data = rbind(df1, df2)
  trees = lidR::las_update(trees)

  return(trees)
}


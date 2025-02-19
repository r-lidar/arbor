#' Extend the trees to the DTM
#'
#' Extend the trees to the DTM. It fits circles using a RANSAC approach to measure the bottom diameter
#' of each tree. In practice, fitting one circle on a low slice of a stem works only in easy contexts with
#' vertical and perfectly segmented wood, and no noise or mistakes. To be more robust, the method fits numerous
#' circles at different heights and with different thicknesses.
#'
#' @param trees LAS object from lidR
#' @param dtm SpatRaster from terra. The DTM.
#' @param circles data.frame produced by \link{measure_diameters}
#' @param extra_height The trees can be prolonged below the ground.
#' @export
tree_extensions = function(trees, dtm, circles, extra_height = 0.15)
{
  pointID <- treeID <- Z <- NULL

  extensions = list()
  trees@data$pointID = 1:lidR::npoints(trees)

  for (i in 1:nrow(circles))
  {
    id  = circles$treeID[i]
    x = circles$center_x[i]
    y = circles$center_y[i]
    z = circles$center_z[i]
    r = circles$radius[i]
    HAG = circles$center_hag[i]
    circle = list(center_x = x, center_y = y, z = z, radius = r)

    loc = matrix(c(x,y), ncol = 2)
    zgnd = terra::extract(dtm, loc, method = "bilinear")
    hcyl = as.numeric(z-zgnd+extra_height)

    cat("Tree", id, "\n")

    tt = lidR::filter_poi(trees, treeID == id, Z < z)
    trees$pointID[tt$pointID] = 0

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

    #x = plot(tt)
    #plot(LAS(extension), add = x)
  }

  extensions = do.call(rbind, extensions)
  extensions$anisotropy = 1
  extensions$wood = TRUE
  data.table::setDT(extensions)

  extensions = lidR::LAS(extensions)

  trees = lidR::filter_poi(trees, pointID > 0)
  trees = weld_extension(trees, extensions)

  return(trees)
}

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


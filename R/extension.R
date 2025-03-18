#' Extend the trees to the DTM
#'
#' Extend the trees to the DTM. It uses fit circle produced by \link{measure_diameters} and remove,
#' for each tree, every point below this circle and replace by a synthetic cylinder. This ensure a
#' proper cleaning of the bottom of the tree, the extension to the ground and a more robust QSM.
#'
#' @param trees LAS object from lidR
#' @param dtm SpatRaster from terra. The DTM.
#' @param circles data.frame produced by \link{measure_diameters}
#' @param extra_height The trees can be prolonged below the ground.
#' @export
tree_extensions = function(trees, dtm, circles, extra_height = 0.15)
{
  pointID <- treeID <- Z <- . <- NULL

  extensions = list()
  trees@data$pointID = 1:lidR::npoints(trees)

  data.table::setindex(trees@data, treeID)

  # Initialize progress bar
  pb <- utils::txtProgressBar(min = 0, max = nrow(circles), style = 3)

  for (i in 1:nrow(circles))
  {
    id  = circles$treeID[i]
    x = circles$center_x[i]
    y = circles$center_y[i]
    z = circles$center_z[i]
    r = circles$radius[i]
    HAG = circles$center_hag[i]
    circle = list(center_x = x, center_y = y, z = z, radius = r)

    loc = matrix(c(x, y), ncol = 2)
    zgnd = terra::extract(dtm, loc, method = "bilinear")
    hcyl = as.numeric(z - zgnd + extra_height)

    ids = trees@data[.(id), on = "treeID"][Z < z]$pointID
    trees$pointID[ids] = 0

    extension = generate_cylinder_points(circle, height = hcyl)
    lidR::quantize(extension$X, trees@header[["X scale factor"]], trees@header[["X offset"]])
    lidR::quantize(extension$Y, trees@header[["Y scale factor"]], trees@header[["Y offset"]])
    lidR::quantize(extension$Z, trees@header[["Z scale factor"]], trees@header[["Z offset"]])
    extension$treeID = id
    extensions[[as.character(id)]] = extension

    # Update progress bar
    utils::setTxtProgressBar(pb, i)
  }

  # Close progress bar
  close(pb)

  extensions = do.call(rbind, extensions)
  extensions$anisotropy = 1
  extensions$wood = TRUE
  data.table::setDT(extensions)

  extensions = lidR::LAS(extensions)

  trees = lidR::filter_poi(trees, pointID > 0)
  trees@data$pointID = NULL
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


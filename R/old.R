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
#' @noRd
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


#' Find a seed for each tree
#'
#' To segment individual tree instances, the first step is to identify seed points for each tree.
#' This is achieved by extracting a slice of points near the ground. The wood/foliage
#' semantic segmentation must be performed first using \link{segment_foliage} to ensure
#' that wood points are correctly identified. The method uses only points that are labelled 'wood'.
#' The method then applies connected component clustering followed by RANSAC circle fitting to reliably
#' assign a single seed to each tree.
#'
#' @param las A LAS object from lidR.
#' @param ... Unused. Additional parameters beyond \code{...} should generally remain unchanged,
#'   except in edge cases.
#' @param max_diameter Maximum expected tree diameter (in meters). Used to filter out invalid
#'   RANSAC-fitted circles.
#' @param slice_seeds_at A numeric vector of two values defining the height range (in meters)
#'   for slicing the point cloud.
#' @param res Resolution for connected component clustering.
# @param smooth Smoothing radius (in meters) applied to the slice to improve the clustering process.
#' @param tree_spacing numeric. If the scene is a plantation, providing the spacing of the trees
#' guarantee much better seeds finding by providing prior information. We now know that trees are
#' regularly spaced.
#' @noRd
find_seeds = function(las, slice_seeds_at = c(0.5, 0.8), ..., max_diameter = 0.6, res = 0.025, tree_spacing = NULL)
{
  treeID <- clusterID <- X <- Y <- Z <-  hag <- hag_max <- hag_min <- foliage <- NULL

  if (!is.null(tree_spacing) && tree_spacing > 0)
  {
    res = tree_spacing/4
  }

  # The point cloud must have hag, anisotropy and foliage computed
  attributes = names(las)
  stopifnot("foliage" %in% attributes)
  stopifnot("hag" %in% attributes)

  cat("Finding seeds with a clusetring approach\n")

  seed = lidR::filter_poi(las, hag > slice_seeds_at[1], hag < slice_seeds_at[2], foliage == FALSE)
  seed = smooth3d(seed, 0.03)
  seed$Z = seed$Z * 0.01
  seed = lidR::connected_components(seed, res, 10)
  seed = lidR::filter_poi(seed, clusterID != 0)
  seed$Z = seed$Z * 100
  #plot(seed, color = "clusterID", pal = pastel.colors(500)) |> add_dtm3d(dtm)

  fit_circle_to_seed = function(id, max_radius)
  {
    cl = seed[seed$clusterID == id]
    if (lidR::npoints(cl) < 10) return(NULL)
    circle = ransac_circle(cl, num_iterations = 400)
    valid  = is.valid.circle(circle$radius, circle$covered_arc_degree, circle$percentage_inlier*100, circle$percentage_inside*100)
    if (valid) return(data.frame(X = circle$center_x, Y = circle$center_y, Z = circle$z, R = circle$radius, id = id))
    else return(NULL)
  }

  cat("Fitting RANSAC circles to each seed\n")

  circles = lapply(unique(seed$clusterID), fit_circle_to_seed, max_radius = max_diameter/2)
  circles = do.call(rbind, circles)

  sfcircles = NULL
  if (!is.null(circles))
  {
    sfcenters = sf::st_as_sf(circles, coords = c("X", "Y", "Z"))
    sfcircles = sf::st_buffer(sfcenters, circles$R*1.20, nQuadSegs = 10)
  }

  f = function(x,y,z)
  {
    Z = c(round(min(z), 3), round(max(z), 3))
    dZ = diff(Z)

    bottom = z < Z[2] - dZ/2
    top = z > Z[1] + dZ/2

    X = c(round(stats::median(x[bottom]), 3), round(stats::median(x[top]), 3))
    Y = c(round(stats::median(y[bottom]), 3), round(stats::median(y[top]), 3))

    return(list(X = X, Y = Y,Z = Z))
  }

  cat("Seed correction with RANSAC circles\n")

  seeds = seed@data[, f(X,Y,Z), by = clusterID]
  seeds = stats::na.omit(seeds)

  sfseeds = sf::st_as_sf(seeds, coords = c("X", "Y", "Z"))

  extraseed = NULL

  if (!is.null(sfcircles))
  {
    intersect = sf::st_intersects(sfcircles, sfseeds)
    ii = lapply(intersect, length)
    ii = which(ii > 2)

    for (iii in ii)
    {
      ids = sfseeds[intersect[[iii]],]$clusterID
      sfseeds$clusterID[sfseeds$clusterID %in% ids] = ids[1]
    }

    cat("Seed generation from RANSAC circles\n")

    sfcenters = sfcenters[sfcenters$id %in% sfseeds$clusterID,]
    extraseed = sf::st_buffer(sfcenters, sfcenters$R*0.9, nQuadSegs = 3)
    extraseed$Z = sf::st_coordinates(sfcenters)[,3]
    extraseed = sf::st_cast(extraseed, "POINT")
    extraseed$R = NULL
    coord = as.data.frame(sf::st_coordinates(extraseed))
    coord$Z = extraseed$Z
    coord$clusterID = extraseed$id
    extraseed = sf::st_as_sf(coord, coords = c("X", "Y", "Z"))
  }

  fullseed = rbind(sfseeds, extraseed)



  #circ = circles[61,]
  #see = sfseeds[intersect[[61]],]
  #plot(sf::st_geometry(circ), axes = T)
  #plot(see, add = T, pch = as.numeric(as.factor(see$clusterID)), col = "black")
  #plot(seed[seed$clusterID %in% see$clusterID], color = "clusterID")

  names(fullseed)[1] = "treeID"

  #plot(seed, color = "clusterID", pal = pastel.colors(500)) |> add_dtm3d(dtm) |> add_treetops3d(sfseeds, radius = 0.05)

  return(fullseed)
}

#' Post production tool
#'
#' Apply this function to fix some issue in the instance segmentation. If you have no issue but apply
#' this function this should not have any effect. This function fits circles on each tree at different
#' levels and tries to figure out if two segmented trees are actually the same one badly segmented because
#' of bad seeds. In some cases tt is possible that one tree have two different seeds leading to two
#' different instance. This function tries to fix that.
#'
#' @param las LAS object from lidR
#' @param max_height maximum height above the lowest point of a tree to fit circles
#' @param slice_thickness slice thickness to fit circles
#' @param max_diameter a threshold to protect against bad circle fitting. If a circle is bigger than
#' that it is not valid. Set it the maximum expected diameter of a tree.
#' @noRd
fix_split_trees = function(las, max_height = 4, slice_thickness = 0.25, max_diameter = 0.6)
{
  treeID <- X <- Y <- Z <-  hag <- hag_max <- hag_min <- foliage <- NULL

  attributes = names(las)
  stopifnot("foliage" %in% attributes)
  stopifnot("hag" %in% attributes)

  max_radius = max_diameter/2

  # Define a function to fit circle
  fit_circle_to_seed = function(id)
  {
    #cat("Tree", id, ":")
    keep = slice$treeID == id
    keep[is.na(keep)] = FALSE
    cl = slice[keep]
    if (lidR::npoints(cl) < 10)
    {
      #cat(" not enought points\n")
      return(NULL)
    }
    circle = ransac_circle(cl, num_iterations = 400, inlier_threshold = 0.02)
    valid = is.valid.circle(circle$radius, circle$covered_arc_degree, circle$percentage_inlier*100, circle$percentage_inside*100)
    if (valid)
    {
      #cat("\n")
      return(data.frame(X = circle$center_x, Y = circle$center_y, Z = circle$z, R = circle$radius, id = id))
    }
    else
    {
      #cat(" invalid\n")
      return(NULL)
    }
  }

  zmin = min(las$hag)
  offsets = seq(0, max_height,slice_thickness)

  # For each slice apply the fit to each tree
  # if two circle intersect they are from the same tree
  for (k in 1:length(offsets))
  {
    slice = lidR::filter_poi(las, hag > zmin+offsets[k], hag < zmin+offsets[k]+slice_thickness, foliage == FALSE)

    #x = plot(slice, color = "treeID")
    circles = lapply(unique(slice$treeID), fit_circle_to_seed)
    circles = do.call(rbind, circles)

    if (is.null(circles)) next

    #for (i in 1:nrow(circles)) add_circle3d(x, circles$X[i], circles$Y[i], circles$R[i], 0)

    circles = sf::st_as_sf(circles, coords = c("X", "Y"))
    circles = sf::st_buffer(circles, circles$R*1.20)

    intersect = sf::st_intersects(circles)
    intersect = Filter(function(x) length(x) > 1, intersect)
    intersect = unique(intersect)

    if (length(intersect) == 0) next

    for (i in intersect)
    {
      intersecting_circles = circles[i,]
      ids = intersecting_circles$id

      intersecting_circles = sf::st_geometry(intersecting_circles)

      if (length(intersecting_circles) > 2)
        warning("More than 2 rings intersecting. This case is not handled yet.")

      area1 = mean(sf::st_area(intersecting_circles))
      area2 = sf::st_area(sf::st_intersection(intersecting_circles[1],intersecting_circles[2]))


      if (area2 > area1*0.5)
      {
        #plot(filter_poi(las, treeID %in% ids), color = "treeID", axis = T)
        #plot(filter_poi(las, treeID %in% ids), color = "foliage")
        cat("Combining trees", ids, "into tree", ids[1], "\n")
        las$treeID[las$treeID %in% ids] = ids[1]
      }
    }
  }

  return(las)
}


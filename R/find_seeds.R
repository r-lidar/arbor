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
#' @export
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
    valid  = is.valid.circle(circle$radius, circle$covered_arc_degree, circle$percentage_inlier*100)
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

seeds_to_sf_point = function(seed)
{

}

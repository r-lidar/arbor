#' Generate a QSM from a single tree point cloud
#'
#' This function processes a tree point cloud to generate a Quantitative Structure Model (QSM) using
#' the AdTree algorithm. The original AdTree software is available at \url{https://github.com/tudelft3d/AdTree}.
#' The version shipped with the package is a modified version with more control on the parameters.
#' The source code is available at \url{https://github.com/Jean-Romain/AdTree}.
#'
#'
#' @param tree A `LAS` object containing a single tree point cloud with required attributes. Only
#' the point labelled as wood will be used for QSM.
#' @param alpha Numeric. Main skeleton simplification with Douglas-Peukey algorithm. Default is 0.5.
#' Was 1 in the original version of AdTree.
#' @param subtree Numeric. Eliminate small edges and keep the main skeleton. Default is 0.02.
#' Was 0.019 in the original version of AdTree.
#' @param min_radius Numeric. The minimum branch radius for the QSM. Reduces the complexity of the QSM
#' by removing very small ramification. Default is 0.003. Does not exist in the original version of AdTree.
#'
#' @details
#' The function extracts non-foliage points, writes them to a temporary file, and then runs the AdTree
#' algorithm to generate a 3D mesh representation of the tree structure. The resulting mesh is loaded
#' and returned as a `Morpho` object using `Morpho::obj2mesh()`
#'
#'
#' @return A `Morpho` mesh object representing the reconstructed tree structure.
#'
#' @examples
#' \dontrun{
#' library(lidR)
#' f <- system.file("extdata", "tree.laz", package="lidRtls")
#' tree <- readTLS(f)
#' qsm <- qsm(tree)
#' plot(tree, color = "foliage", pal = c("chocolate4", "darkgreen"), clear_artifacts = FALSE)
#' rgl::shade3d(qsm, col="brown4")
#' }
#' @export
qsm = function(tree, alpha = 0.8, subtree = 0.02, min_radius = 0.003)
{
  . <- X <- Y <- Z <- foliage <- NULL

  attributes = names(tree)
  stopifnot("anisotropy" %in% attributes)
  stopifnot("foliage" %in% attributes)
  stopifnot("treeID" %in% attributes)

  if (length(unique(tree$treeID)) != 1) stop("The point cloud must contain a single tree", call. = FALSE)

  ofile = tempfile(fileext = ".xyz")
  odir = tempdir()
  iobj = tools::file_path_sans_ext(ofile)
  iobj = paste0(iobj, "_branches.obj")

  no_foliage = lidR::filter_poi(tree, foliage == FALSE)
  xyz = no_foliage@data[, .(X,Y,Z)]
  data.table::fwrite(xyz, ofile, sep = " ", col.names = FALSE)

  os = tolower(Sys.info()["sysname"])
  bin = paste0("AdTree-", os)

  if (os == "darwin")
    stop("A version of AdTree is avaiable for Mac but because I can't test it I did not included it. Please contact us.")

  adtree = paste0(system.file("bin", "", package = "lidRtls"), "/", bin)

  if (os == "linux")
    adtree <- gsub(" ", "\\\\ ", adtree)

  if (os == "windows")
    adtree = paste0('"', adtree, '.exe"')

  args = paste0("-radius ", min_radius, " -alpha ", alpha, " -subtree ", subtree)
  cmd = paste(adtree, ofile, odir, args)

  suppressWarnings(system(cmd, TRUE))

  if (!file.exists(iobj)) stop("Failed to produced QSM .obj file")

  mesh = Morpho::obj2mesh(iobj)
  mesh
}


#' Plot and Add 3D QSM Objects
#'
#' These functions allow you to visualize 3D Quantitative Structure Models (QSM) using `rgl`.
#' `plot_qsm3d()` initializes a new 3D plot and displays the QSM, while `add_qsm3d()` adds another QSM to the existing plot.
#'
#' @param qsm A mesh object representing the QSM, produced by \link{qsm}.
#' @param bottom_to_zero Logical. If `TRUE`, shifts the QSM so that the lowest z-coordinate is at zero. Default is `TRUE`.
#' @param x A numeric vector of length 2 containing x and y offsets to align QSMs and point cloud.
#'
#' @return `plot_qsm3d()` returns an invisible numeric vector of length 2 containing the original minimum x and y coordinates.
#' `add_qsm3d()` does not return a value.
#'
#' @md
#' @export
plot_qsm3d = function(qsm, bottom_to_zero = TRUE)
{
  x = min(qsm$vb[1,])
  y = min(qsm$vb[2,])
  z = min(qsm$vb[3,])
  qsm$vb[1,] = qsm$vb[1,] - x
  qsm$vb[2,] = qsm$vb[2,] - y
  if (bottom_to_zero) qsm$vb[3,] = qsm$vb[3,] - z
  rgl::open3d()
  rgl::shade3d(qsm, col="brown4")
  return(invisible(c(x,y)))
}

#' @rdname plot_qsm3d
#' @export
add_qsm3d = function(x, qsm)
{
  qsm$vb[1,] =   qsm$vb[1,] - x[1]
  qsm$vb[2,] =   qsm$vb[2,] - x[2]
  rgl::shade3d(qsm, col="brown4")
}

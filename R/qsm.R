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
qsm = function(tree, alpha = 0.5, subtree = 0.02, min_radius = 0.003)
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

#' @method plot qsf
#' @export
#' @rdname plot
plot.qsf = function(x, color = "cyl_ID", pal = c("blue", "green", "yellow", "orange", "red"))
{
  meshes = lapply(x, as_mesh, color, pal)
  mesches = Filter(Negate(is.null), meshes)
  rgl::open3d()
  rgl::bg3d("black")
  rgl::shapelist3d(mesches)
}


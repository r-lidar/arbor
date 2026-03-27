#' @method plot qsf
#' @export
#' @rdname plot
plot.qsf = function(x, ..., color = "cyl_ID", pal = c("blue", "green", "yellow", "orange", "red"), add = NULL)
{
  if (!is.null(add)) {
    tx = add[1]; ty = add[2]; tz = 0
  } else {
    tx = min(sapply(x, function(qsm) min(qsm$startX)))
    ty = min(sapply(x, function(qsm) min(qsm$startY)))
    tz = min(sapply(x, function(qsm) min(qsm$startZ)))
  }

  x = lapply(x, function(qsm)
  {
    qsm$startX <- qsm$startX - tx
    qsm$startY <- qsm$startY - ty
    qsm$startZ <- qsm$startZ - tz
    qsm$endX   <- qsm$endX - tx
    qsm$endY   <- qsm$endY - ty
    qsm$endZ   <- qsm$endZ - tz
    qsm
  })

  meshes = lapply(x, as_mesh, color, pal)
  mesches = Filter(Negate(is.null), meshes)
  rgl::open3d()
  rgl::bg3d("black")
  rgl::shapelist3d(mesches)
  lidR:::.pan3d(2)
}


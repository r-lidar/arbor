# @file qsf_plot.R
# Project: Arbor
#
# Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' @method plot qsf
#' @export
#' @rdname plot
plot.qsf = function(x, ..., color = "branch_order", pal = "auto", add = NULL)
{
  default_pal = c("blue", "green", "yellow", "orange", "red")
  branch_order_pal <- c("#552203","#ad3a01","#e59b16","#fdd63b","#9acd56", "#59ad4b","#1e490e")
  quality_pal = c("blue", "green", "yellow", "orange", "red")

  if (identical(pal, "auto")) {
    pal <- switch(color, branch_order = branch_order_pal, quality = quality_pal, default_pal)
  }

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
    #qsm$startZ <- qsm$startZ - tz
    qsm$endX   <- qsm$endX - tx
    qsm$endY   <- qsm$endY - ty
    #qsm$endZ   <- qsm$endZ - tz
    qsm
  })

  meshes = lapply(x, as_mesh, color, pal)
  mesches = Filter(Negate(is.null), meshes)
  if (is.null(add)) rgl::open3d()
  rgl::bg3d("black")
  rgl::shapelist3d(mesches)
  lidR:::.pan3d(2)
  return(c(tx, ty))
}


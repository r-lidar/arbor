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
plot.qsf <- function(x, ..., color = "branch_order", pal = "auto", add = NULL) {
  
  # Delegate everything to as_mesh.qsf
  meshes <- as_mesh.qsf(x, offset = add, color = color, pal = pal)
  applied_offset <- attr(meshes, "offset")

  if (is.null(add)) rgl::open3d()
  
  rgl::bg3d("black")
  rgl::shapelist3d(meshes)
  lidR:::.pan3d(2)
  
  return(invisible(applied_offset))
}


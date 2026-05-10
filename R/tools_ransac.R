# @file tools_ransac.R
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

ransac_circle <- function(points, num_iterations = 100, inlier_threshold = 0.01)
{
  if (methods::is(points, "LAS")) points = sf::st_coordinates(points)
  stopifnot(is.matrix(points), ncol(points) == 3L, nrow(points) > 3L)
  return(ransac_circle_cpp(points, num_iterations, inlier_threshold))
}

#' Add 3D Circles to an rgl Plot
#'
#' This function adds multiple 3D circles to an existing rgl plot.
#'
#' @param x A numeric vector of length 2, representing the reference point (x, y).
#' @param center_x A numeric vector of x-coordinates for the circle centers.
#' @param center_y A numeric vector of y-coordinates for the circle centers.
#' @param radius A numeric vector of radii for the circles.
#' @param height A numeric vector specifying the z-coordinate (height) for each circle.
#'
#' @return No return value. This function modifies the rgl plot by adding 3D circles.
#' @noRd
add_circles3d <- function(x, center_x, center_y, radius, height)
{
  theta <- seq(0, 2 * pi, length.out = 50)
  cos_theta <- cos(theta)
  sin_theta <- sin(theta)

  for (i in seq_along(center_x))
  {
    xc <- center_x[i]
    yc <- center_y[i]
    r <- radius[i]
    h <- height[i]

    xx <- xc - x[1] + r * cos_theta
    yy <- yc - x[2] + r * sin_theta
    zz <- rep(h, length(theta))

    rgl::lines3d(xx, yy, zz, col = "red", lwd = 3)
  }
}


add_circles3d <- function(x, center_x, center_y, radius, height)
{
  theta <- seq(0, 2 * pi, length.out = 50)
  cos_theta <- cos(theta)
  sin_theta <- sin(theta)

  rings <- list()

  for (i in seq_along(center_x))
  {
    xc <- center_x[i]
    yc <- center_y[i]
    r <- radius[i]
    h <- height[i]

    xx <- center_x - x[1] + radius * cos_theta
    yy <- center_y - x[2] + radius * sin_theta
    zz <- rep(height, length(theta))

    rings[[i]] <- rgl::lines3d(xx, yy, zz)
  }

  rgl::shade3d(rgl::shapelist3d(cyls, plot = FALSE), col = "red")
}

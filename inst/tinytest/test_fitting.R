# Test case 7: RANSAC on noisy circle points
fit = arbor:::fit_circloid_cpp

set.seed(123)
n = 500
theta <- seq(0, 2 * pi, length.out = n)
xc = 12
yc = 18.5
r = 2
circle_points <- matrix(c(xc + r * cos(theta) + rnorm(n, 0, 0.1),
                          yc + r * sin(theta) + rnorm(n, 0, 0.1),
                          rep(0, n)), ncol = 3, byrow = FALSE)


a <- 4.0
b <- 1.5

ellipse_points <- matrix(c(
  xc + a * cos(theta) + rnorm(n, 0, 0.1),   # X coordinates
  yc + b * sin(theta) + rnorm(n, 0, 0.1), # Y coordinates
  rep(0, n)                                 # Z coordinates (remains 0)
), ncol = 3, byrow = FALSE)


r = 3+0.5*cos(theta)+sin(2*theta) + runif(n, -0.25, 0.25)
x = xc + r*cos(theta)
y = yc + r*sin(theta)
z = 0
circloid_points = cbind(x,y,z)


# ================
# Check circle
# ================

res <- fit(circle_points)

if (FALSE)
{
  plot(circle_points, asp = 1, main = res$shape_type)
  lines(res$nodes, lwd = 2, col = "red")
  points(res$center_x, res$center_y, pch = 3, cex = 2)
}

expect_equal(res$center_x, xc, tolerance = 0.005)
expect_equal(res$center_y, yc, tolerance = 0.005)
expect_equal(res$radius, r, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 360)


# ================
# Check circle
# ================

res <- fit(ellipse_points)

if (FALSE)
{
  plot(ellipse_points, asp = 1, main = res$shape_type)
  lines(res$nodes, lwd = 2, col = "red")
  points(res$center_x, res$center_y, pch = 3, cex = 2)
}

expect_equal(res$center_x, xc, tolerance = 0.005)
expect_equal(res$center_y, yc, tolerance = 0.005)
expect_equal(res$radius, 2.51, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 360)

# ================
# Check circloid
# ================

res <- fit(circloid_points, tolerance = 0.05)

if (FALSE)
{
  plot(circloid_points, asp = 1, main = res$shape_type)
  lines(res$nodes, lwd = 2, col = "red")
  points(res$center_x, res$center_y, pch = 3, cex = 2)
}
expect_equal(res$center_x, xc, tolerance = 0.05)
expect_equal(res$center_y, yc, tolerance = 0.05)
expect_equal(res$radius, 3.1, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 360)


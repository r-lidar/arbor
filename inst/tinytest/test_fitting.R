# Test case 7: RANSAC on noisy circle points
fit = arbor:::fit_circloid_cpp

set.seed(123)
n = 500
theta <- seq(0, 2 * pi, length.out = n)
xc = 12
yc = 18.5
r = 2
a <- 4.0
b <- 1.5

circle_points <- matrix(c(xc + r * cos(theta) + rnorm(n, 0, 0.1),
                          yc + r * sin(theta) + rnorm(n, 0, 0.1),
                          rep(0, n)), ncol = 3, byrow = FALSE)


ellipse_points <- matrix(c(xc + a * cos(theta) + rnorm(n, 0, 0.1),
                           yc + b * sin(theta) + rnorm(n, 0, 0.1),
                           rep(0, n)), ncol = 3, byrow = FALSE)

theta = seq(0, pi, length.out = n)

hcircle_points <- matrix(c(xc + r * cos(theta) + rnorm(n, 0, 0.1),
                           yc + r * sin(theta) + rnorm(n, 0, 0.1),
                           rep(0, n)), ncol = 3, byrow = FALSE)


hellipse_points <- matrix(c(xc + a * cos(theta) + rnorm(n, 0, 0.1),
                            yc + b * sin(theta) + rnorm(n, 0, 0.1),
                            rep(0, n)), ncol = 3, byrow = FALSE)

theta <- seq(0, 2 * pi, length.out = n)
r = 3+0.5*cos(theta)+sin(2*theta) + runif(n, -0.25, 0.25)
x = xc + r*cos(theta)
y = yc + r*sin(theta)
z = 0
circloid_points = cbind(x,y,z )

theta = seq(0, pi, length.out = n)
r = 3+0.5*cos(theta)+sin(2*theta) + runif(n, -0.25, 0.25)
x = xc + r*cos(theta)
y = yc + r*sin(theta)
z = 0
hcircloid_points = cbind(x,y,z )



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
# Check half circle
# ================

res <- fit(hcircle_points, tolerance = 0.05)

if (FALSE)
{
  plot(hcircle_points, asp = 1, main = res$shape_type)
  lines(res$nodes, lwd = 2, col = "red")
  points(res$center_x, res$center_y, pch = 3, cex = 2)
}

expect_equal(res$center_x, xc, tolerance = 0.005)
expect_equal(res$center_y, yc, tolerance = 0.005)
expect_equal(res$radius, r, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 180)


# ================
# Check ellipse
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
# Check half ellipse
# ================

res <- fit(hellipse_points)

if (FALSE)
{
  plot(hellipse_points, asp = 1, main = res$shape_type)
  lines(res$nodes, lwd = 2, col = "red")
  points(res$center_x, res$center_y, pch = 3, cex = 2)
}

expect_equal(res$center_x, xc, tolerance = 0.005)
expect_equal(res$center_y, yc, tolerance = 0.005)
expect_equal(res$radius, 2.51, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 180)

# ================
# Check half ellipse
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


# ================
# Check half circloid
# ================

res <- fit(hcircloid_points, tolerance = 0.05)

if (FALSE)
{
  plot(hcircloid_points, asp = 1, main = res$shape_type)
  lines(res$nodes, lwd = 2, col = "red")
  points(res$center_x, res$center_y, pch = 3, cex = 2)
}



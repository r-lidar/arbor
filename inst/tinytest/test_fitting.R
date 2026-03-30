<<<<<<< HEAD
fit = arbor:::fit_circloid_cpp

show = function(pt, res)
{
  inliers = pt[res$inliers,]
  plot(pt, asp = 1, main = res$shape_type)
  points(res$center_x, res$center_y, pch = 3, cex = 2)
  points(inliers, col = "blue", pch = 18)
  lines(res$nodes, lwd = 2, col = "red")
  symbols(res$center_x, res$center_y, circles = res$radius, inches = FALSE, add = T, fg = "purple")
  symbols(res$center_x, res$center_y, circles = res$radius+tol, inches = FALSE, add = T, fg = "purple")
  symbols(res$center_x, res$center_y, circles = res$radius-tol, inches = FALSE, add = T, fg = "purple")
}

show3d = function(pt, res)
{
  inliers = pt[res$inliers,]
  rgl::points3d(pt)
  rgl::points3d(res$center_x, res$center_y, res$center_z, size = 6, col = "red")
  rgl::points3d(inliers, col = "blue", pch = 18)
  rgl::lines3d(res$nodes, lwd = 2, col = "red")
}


=======
# Test case 7: RANSAC on noisy circle points
fit = arbor:::fit_circloid_cpp

>>>>>>> main
set.seed(123)
n = 500
theta <- seq(0, 2 * pi, length.out = n)
xc = 12
yc = 18.5
r = 2
<<<<<<< HEAD
a <- 4.0
b <- 1.5

theta_rot <- pi / 4  # 45 degrees
Rx <- matrix(c(
  1, 0, 0,
  0, cos(theta_rot), -sin(theta_rot),
  0, sin(theta_rot),  cos(theta_rot)
), ncol = 3, byrow = TRUE)

circle_points <- matrix(c(xc + r * cos(theta) + rnorm(n, 0, 0.1),
                          yc + r * sin(theta) + rnorm(n, 0, 0.1),
                          runif(n, 0, 2)), ncol = 3, byrow = FALSE)


ellipse_points <- matrix(c(xc + a * cos(theta) + rnorm(n, 0, 0.1),
                           yc + b * sin(theta) + rnorm(n, 0, 0.1),
                           runif(n, 0, 2)), ncol = 3, byrow = FALSE)

theta = seq(0, pi, length.out = n)

hcircle_points <- matrix(c(xc + r * cos(theta) + rnorm(n, 0, 0.1),
                           yc + r * sin(theta) + rnorm(n, 0, 0.1),
                           runif(n, 0, 2)), ncol = 3, byrow = FALSE)


hellipse_points <- matrix(c(xc + a * cos(theta) + rnorm(n, 0, 0.1),
                            yc + b * sin(theta) + rnorm(n, 0, 0.1),
                            runif(n, 0, 2)), ncol = 3, byrow = FALSE)

rcircle_points <- circle_points %*% t(Rx)


theta <- seq(0, 2 * pi, length.out = n)
R = 3+0.5*cos(theta)+sin(2*theta) + runif(n, -0.25, 0.25)
x = xc + R*cos(theta)
y = yc + R*sin(theta)
z = runif(n, 0, 2)
circloid_points = cbind(x,y,z )

rcircloid_points <- circloid_points %*% t(Rx)

theta = seq(0, pi, length.out = n)
R = 3+0.5*cos(theta)+sin(2*theta) + runif(n, -0.25, 0.25)
x = xc + R*cos(theta)
y = yc + R*sin(theta)
z = 0
hcircloid_points = cbind(x,y,z )

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

pt = circle_points
tol = 0.07
res <- fit(pt, tolerance = tol)

if (FALSE) show(pt, res)

expect_equal(res$center_x, xc, tolerance = 0.005)
expect_equal(res$center_y, yc, tolerance = 0.005)
expect_equal(res$radius, r, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 360)


pt = circle_points
tol = 0.15
res <- fit(pt, tolerance = tol)

if (FALSE) show(pt, res)

expect_equal(res$center_x, xc, tolerance = 0.005)
expect_equal(res$center_y, yc, tolerance = 0.005)
expect_equal(res$radius, r, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 360)

# ================
# Check half circle
# ================

pt = hcircle_points
tol = 0.07
res <- fit(pt, tolerance = tol)

if (FALSE) show(pt, res)

expect_equal(res$center_x, xc, tolerance = 0.005)
expect_equal(res$center_y, yc, tolerance = 0.005)
expect_equal(res$radius, r, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 180)

pt = hcircle_points
tol = 0.15
res <- fit(pt, tolerance = tol)

if (FALSE) show(pt, res)

expect_equal(res$center_x, xc, tolerance = 0.005)
expect_equal(res$center_y, yc, tolerance = 0.005)
expect_equal(res$radius, r, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 180)

# ================
# Check rotated circle
# ================

pt = rcircle_points
tol = 0.15
res <- fit(pt, tolerance = tol, from = c(0,0,0), to = c(0, -1, 1))

if (FALSE) show3d(pt, res)

expect_equal(res$center_x, xc, tolerance = 0.005)
expect_equal(res$center_y, 12.38, tolerance = 0.005)
expect_equal(res$radius, r, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 360)

# ================
# Check ellipse
# ================

pt = ellipse_points
tol = 0.08
res <- fit(pt, tolerance = tol)

if (FALSE) show(pt, res)

expect_equal(res$center_x, xc, tolerance = 0.005)
expect_equal(res$center_y, yc, tolerance = 0.005)
expect_equal(res$radius, 2.41, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 360)


pt = ellipse_points
tol = 0.15
res <- fit(pt, tolerance = tol)

if (FALSE) show(pt, res)

expect_equal(res$center_x, xc, tolerance = 0.005)
expect_equal(res$center_y, yc, tolerance = 0.005)
expect_equal(res$radius, 2.41, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 360)

# ================
# Check half ellipse
# ================

pt = hellipse_points
tol = 0.08
res <- fit(pt, tolerance = tol)

if (FALSE) show(pt, res)

expect_equal(res$center_x, xc, tolerance = 0.005)
expect_equal(res$center_y, yc, tolerance = 0.005)
expect_equal(res$radius, 2.41, tolerance = 0.05)
expect_equal(res$covered_arc_degree, 180)


# ================
# Check circloid
# ================

pt = circloid_points
tol = 0.07
res <- fit(pt, tolerance = tol)

if (FALSE) show(pt, res)

expect_equal(res$center_x, xc, tolerance = 0.05)
expect_equal(res$center_y, yc, tolerance = 0.05)
expect_equal(res$radius, 3.1, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 360)


pt = circloid_points
tol = 0.15
res <- fit(pt, tolerance = tol)
if (FALSE) show(pt, res)

expect_equal(res$center_x, xc, tolerance = 0.05)
expect_equal(res$center_y, yc, tolerance = 0.05)
expect_equal(res$radius, 3.1, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 360)

# ================
# Check half circloid
# ================

pt = hcircloid_points
tol = 0.2
res <- fit(pt, tolerance = tol)

if (FALSE)


expect_equal(res$center_x, 13.4, tolerance = 0.01)
expect_equal(res$center_y, 19.6, tolerance = 0.01)
expect_equal(res$radius, 2.43, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 220)

# ================
# Check rotated circloid
# ================

pt = rcircloid_points
tol = 0.15
res <- fit(pt, tolerance = tol, from = c(0,0,0), to = c(0, -1, 1))

if (FALSE) show3d(pt, res)

expect_equal(res$center_x, xc, tolerance = 0.05)
expect_equal(res$center_y, 12.48, tolerance = 0.05)
expect_equal(res$radius, 3.1, tolerance = 0.025)
expect_equal(res$covered_arc_degree, 360)


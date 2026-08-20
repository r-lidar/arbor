fit <- arbor:::fit_circloid_cpp

show = function(pt, res, tol = 0.03)
{
  inliers = pt[res$inliers,]

  plot(pt, asp = 1, main = res$shape_type)
  points(res$center_x, res$center_y, pch = 3, cex = 2)
  points(inliers, col = "blue", pch = 18)
  lines(res$nodes, lwd = 2, col = "red")

  symbols(res$center_x, res$center_y,
          circles = res$radius,
          inches = FALSE, add = TRUE, fg = "purple")

  symbols(res$center_x, res$center_y,
          circles = res$radius + tol,
          inches = FALSE, add = TRUE, fg = "purple")

  symbols(res$center_x, res$center_y,
          circles = res$radius - tol,
          inches = FALSE, add = TRUE, fg = "purple")

  # Draw 36 sectors (10-degree spacing)
  theta <- seq(0, 350, by = 10) * pi / 180

  segments(
    x0 = res$center_x,
    y0 = res$center_y,
    x1 = res$center_x + (res$radius + tol) * cos(theta),
    y1 = res$center_y + (res$radius + tol) * sin(theta),
    col = "grey80"
  )
}


show3d = function(pt, res)
{
  inliers = pt[res$inliers,]

  rgl::points3d(pt)
  rgl::points3d(res$center_x, res$center_y, res$center_z,
                size = 6, col = "red")
  rgl::points3d(inliers, col = "blue", pch = 18)
  rgl::lines3d(res$nodes, lwd = 2, col = "red")
}

set.seed(123)
n <- 500
theta <- seq(0, 2 * pi, length.out = n)
xc <- 12
yc <- 18.5
r <- 2
a <- 4.0
b <- 1.5

theta_rot <- pi / 4
Rx <- matrix(c(
  1, 0, 0,
  0, cos(theta_rot), -sin(theta_rot),
  0, sin(theta_rot),  cos(theta_rot)
), ncol = 3, byrow = TRUE)

set.seed(123)
ex <- rnorm(n, 0, 0.1)
set.seed(321)
ey <- rnorm(n, 0, 0.1)

circle_points <- matrix(c(xc + r * cos(theta) + ex,
                          yc + r * sin(theta) + ey,
                          runif(n, 0, 2)), ncol = 3, byrow = FALSE)

ellipse_points <- matrix(c(xc + a * cos(theta) + ex,
                           yc + b * sin(theta) + ey,
                           runif(n, 0, 2)), ncol = 3, byrow = FALSE)

theta_half <- seq(0, pi, length.out = n)

hcircle_points <- matrix(c(xc + r * cos(theta_half) + ex,
                           yc + r * sin(theta_half) + ey,
                           runif(n, 0, 2)), ncol = 3, byrow = FALSE)

hellipse_points <- matrix(c(xc + a * cos(theta_half) + ex,
                            yc + b * sin(theta_half) + ey,
                            runif(n, 0, 2)), ncol = 3, byrow = FALSE)

rcircle_points <- circle_points %*% t(Rx)

set.seed(42)
R_vals <- 3 + 0.5 * cos(theta) + sin(2 * theta) + runif(n, -0.25, 0.25)
circloid_points <- cbind(xc + R_vals * cos(theta),
                         yc + R_vals * sin(theta),
                         runif(n, 0, 2))

rcircloid_points <- circloid_points %*% t(Rx)

theta_hc <- seq(0, pi, length.out = n)
set.seed(314)
R_hc <- 3 + 0.5 * cos(theta_hc) + sin(2 * theta_hc) + runif(n, -0.25, 0.25)
hcircloid_points <- cbind(xc + R_hc * cos(theta_hc),
                          yc + R_hc * sin(theta_hc),
                          rep(0, n))

disp = FALSE

test_that("fitting detects full circle (tol = 0.07)", {
  res <- fit(circle_points, tolerance = 0.07)

  expect_equal(res$center_x, xc, tolerance = 0.005)
  expect_equal(res$center_y, yc, tolerance = 0.005)
  expect_equal(res$radius, r, tolerance = 0.025)
  expect_equal(res$covered_arc_degree, 360)
  expect_equal(res$shape_type, "circle")


  if (disp) show(circle_points, res)
})

test_that("fitting detects full circle (tol = 0.15)", {
  res <- fit(circle_points, tolerance = 0.15)

  expect_equal(res$center_x, xc, tolerance = 0.005)
  expect_equal(res$center_y, yc, tolerance = 0.005)
  expect_equal(res$radius, r, tolerance = 0.025)
  expect_equal(res$covered_arc_degree, 360)
  expect_equal(res$shape_type, "circle")

  if (disp) show(circle_points, res)
})

test_that("fitting detects half circle (tol = 0.06)", {
  res <- fit(hcircle_points, tolerance = 0.06)

  expect_equal(res$center_x, xc, tolerance = 0.005)
  expect_equal(res$center_y, yc, tolerance = 0.005)
  expect_equal(res$radius, r, tolerance = 0.025)
  expect_equal(res$covered_arc_degree, 200)
  expect_equal(res$shape_type, "circle")

  if (disp) show(hcircle_points, res)
})

test_that("fitting detects half circle (tol = 0.1)", {
  res <- fit(hcircle_points, tolerance = 0.1)

  expect_equal(res$center_x, xc, tolerance = 0.005)
  expect_equal(res$center_y, yc, tolerance = 0.005)
  expect_equal(res$radius, r, tolerance = 0.025)
  expect_equal(res$covered_arc_degree, 200)
  expect_equal(res$shape_type, "circle")

  if (disp) show(hcircle_points, res)
})

test_that("fitting detects rotated circle", {
  res <- fit(rcircle_points, tolerance = 0.15, from = c(0, 0, 0), to = c(0, -1, 1))

  expect_equal(res$center_x, xc, tolerance = 0.005)
  expect_equal(res$center_y, 12.38, tolerance = 0.005)
  expect_equal(res$radius, r, tolerance = 0.025)
  expect_equal(res$covered_arc_degree, 360)
  expect_equal(res$shape_type, "circle")

  if (disp) show3d(rcircle_points, res)
})

test_that("fitting detects full ellipse (tol = 0.08)", {
  res <- fit(ellipse_points, tolerance = 0.08)

  expect_equal(res$center_x, xc, tolerance = 0.005)
  expect_equal(res$center_y, yc, tolerance = 0.005)
  expect_equal(res$radius, 2.41, tolerance = 0.025)
  expect_equal(res$covered_arc_degree, 350)
  expect_equal(res$shape_type, "ellipse")

  if (disp) show(ellipse_points, res)
})

test_that("fitting detects full ellipse (tol = 0.15)", {
  res <- fit(ellipse_points, tolerance = 0.15)

  expect_equal(res$center_x, xc, tolerance = 0.005)
  expect_equal(res$center_y, yc, tolerance = 0.005)
  expect_equal(res$radius, 2.41, tolerance = 0.025)
  expect_equal(res$covered_arc_degree, 360)
  expect_equal(res$shape_type, "ellipse")

  if (disp) show(ellipse_points, res)
})

test_that("fitting detects half ellipse", {
  res <- fit(hellipse_points, tolerance = 0.08)

  expect_equal(res$center_x, xc, tolerance = 0.005)
  expect_equal(res$center_y, 18.75, tolerance = 0.005)
  expect_equal(res$radius, 2.17, tolerance = 0.05)
  expect_true(res$covered_arc_degree >= 180)
  expect_equal(res$shape_type, "ellipse")

  if (disp) show(hellipse_points, res)
})

test_that("fitting detects full circloid (tol = 0.07)", {
  res <- fit(circloid_points, tolerance = 0.07)

  expect_equal(res$center_x, xc, tolerance = 0.05)
  expect_equal(res$center_y, yc, tolerance = 0.05)
  expect_equal(res$radius, 3.1, tolerance = 0.025)
  expect_equal(res$covered_arc_degree, 360)
  expect_equal(res$shape_type, "fourier5")

  if (disp) show(circloid_points, res)
})

test_that("fitting detects full circloid (tol = 0.15)", {
  res <- fit(circloid_points, tolerance = 0.15)

  expect_equal(res$center_x, xc, tolerance = 0.05)
  expect_equal(res$center_y, yc, tolerance = 0.05)
  expect_equal(res$radius, 3.1, tolerance = 0.025)
  expect_equal(res$covered_arc_degree, 360)
  expect_equal(res$shape_type, "fourier5")

  if (disp) show(circloid_points, res)
})

test_that("fitting detects half circloid", {
  res <- fit(hcircloid_points, tolerance = 0.2)

  expect_equal(res$center_x, 13.4, tolerance = 0.02)
  expect_equal(res$center_y, 19.6, tolerance = 0.02)
  expect_equal(res$radius, 2.40, tolerance = 0.05)
  expect_true(res$covered_arc_degree >= 220)
  expect_equal(res$shape_type, "ellipse")

  if (disp) show(hcircloid_points, res)
})

test_that("fitting detects rotated circloid", {
  res <- fit(rcircloid_points, tolerance = 0.15, from = c(0, 0, 0), to = c(0, -1, 1))

  expect_equal(res$center_x, xc, tolerance = 0.05)
  expect_equal(res$center_y, 12.48, tolerance = 0.05)
  expect_equal(res$radius, 3.1, tolerance = 0.025)
  expect_equal(res$covered_arc_degree, 360)
  expect_equal(res$shape_type, "fourier5")


  if (disp) show3d(rcircloid_points, res)
})

test_that("fitting works on real case point cloud", {
  pt <- structure(c(312335.542, 312335.582, 312335.491, 312335.497, 312335.493,
              312335.455, 312335.492, 312335.451, 312335.463, 312335.482, 312335.458,
              312335.477, 312335.436, 312335.413, 312335.487, 312335.469, 312335.488,
              312335.434, 312335.417, 312335.439, 312335.482, 312335.435, 312335.427,
              312335.451, 312335.438, 312335.49, 312335.446, 312335.435, 312335.47,
              312335.497, 312335.495, 312335.453, 312335.473, 312335.465, 312335.469,
              312335.444, 312335.457, 312335.444, 312335.488, 312335.436, 312335.574,
              312335.598, 312335.589, 312335.504, 312335.549, 312335.56, 312335.532,
              312335.568, 312335.521, 312335.514, 312335.565, 312335.585, 312335.546,
              312335.591, 312335.592, 312335.572, 312335.541, 312335.503, 312335.532,
              312335.507, 312335.521, 312335.509, 312335.531, 312335.642, 312335.643,
              312335.618, 312335.694, 312335.622, 312335.693, 312335.652, 312335.634,
              312335.648, 312335.64, 312335.604, 312335.651, 312335.635, 312335.615,
              312335.605, 312335.676, 312335.662, 312335.695, 312335.69, 312335.65,
              312335.616, 312335.618, 312335.677, 312335.698, 312335.62, 312335.671,
              312335.719, 312335.765, 312335.707, 312335.733, 312335.716, 312335.709,
              312335.702, 312335.712, 312335.745, 312335.707, 312335.714, 312335.395,
              312335.398, 312335.399, 312335.411, 312335.416, 312335.428, 312335.436,
              312335.416, 312335.402, 312335.429, 312335.422, 312335.432, 312335.419,
              312335.43, 312335.402, 312335.679, 312335.747, 312335.732, 312335.737,
              312335.712, 312335.723, 312335.74, 312335.737, 312335.752, 312335.771,
              312335.738, 312335.732, 312335.388, 312335.398, 312335.386, 312335.403,
              312335.413, 312335.408, 312335.413, 312335.405, 312335.792, 312335.777,
              312335.757, 312335.761, 312335.757, 312335.77, 312335.772, 312335.753,
              312335.731, 312335.748, 312335.443, 312335.479, 312335.477, 312335.448,
              312335.462, 312335.438, 312335.449, 312335.435, 312335.42, 312335.476,
              312335.443, 312335.482, 312335.567, 312335.669, 312335.688, 312335.675,
              312335.658, 312335.706, 312335.719, 312335.732, 312335.733, 312335.488,
              312335.591, 312335.573, 312335.584, 312335.545, 312335.563, 312335.552,
              312335.524, 312335.681, 312335.636, 312335.662, 312335.68, 312335.383,
              312335.393, 312335.399, 312335.405, 312335.412, 312335.374, 312335.4,
              312335.394, 312335.399, 312335.394, 312335.401, 312335.407, 312335.428,
              312335.402, 312335.427, 312335.423, 312335.461, 312335.424, 312335.423,
              312335.426, 312335.428, 312335.462, 5096506.898, 5096506.892,
              5096506.971, 5096506.976, 5096506.942, 5096506.978, 5096506.945,
              5096506.975, 5096506.946, 5096506.944, 5096506.956, 5096506.939,
              5096506.97, 5096506.995, 5096506.935, 5096506.944, 5096506.945,
              5096506.964, 5096506.981, 5096506.993, 5096506.962, 5096506.949,
              5096506.993, 5096506.978, 5096506.979, 5096506.949, 5096506.971,
              5096506.985, 5096506.946, 5096506.942, 5096506.945, 5096506.982,
              5096506.97, 5096506.976, 5096506.97, 5096506.979, 5096506.977,
              5096506.996, 5096506.952, 5096506.997, 5096506.922, 5096506.938,
              5096506.932, 5096506.936, 5096506.92, 5096506.923, 5096506.918,
              5096506.917, 5096506.923, 5096506.919, 5096506.914, 5096506.921,
              5096506.91, 5096506.944, 5096506.912, 5096506.901, 5096506.938,
              5096506.939, 5096506.955, 5096506.94, 5096506.928, 5096506.936,
              5096506.946, 5096506.964, 5096506.967, 5096506.94, 5096506.997,
              5096506.968, 5096506.971, 5096506.954, 5096506.936, 5096506.941,
              5096506.952, 5096506.942, 5096506.956, 5096506.953, 5096506.947,
              5096506.938, 5096506.943, 5096506.946, 5096506.983, 5096506.998,
              5096506.955, 5096506.931, 5096506.941, 5096506.97, 5096506.963,
              5096506.932, 5096506.954, 5096506.999, 5096506.998, 5096506.972,
              5096506.956, 5096506.976, 5096506.998, 5096506.955, 5096506.972,
              5096506.996, 5096506.971, 5096506.977, 5096507.084, 5096507.053,
              5096507.077, 5096507.006, 5096507.031, 5096507.03, 5096507.041,
              5096507.011, 5096507.06, 5096507.02, 5096507.022, 5096507.006,
              5096507.021, 5096507.008, 5096507.089, 5096507.011, 5096507.038,
              5096507.08, 5096507.056, 5096507.007, 5096507.017, 5096507.039,
              5096507.062, 5096507.092, 5096507.065, 5096507.002, 5096507.033,
              5096507.1, 5096507.159, 5096507.169, 5096507.157, 5096507.149,
              5096507.104, 5096507.192, 5096507.164, 5096507.12, 5096507.16,
              5096507.149, 5096507.167, 5096507.148, 5096507.151, 5096507.157,
              5096507.19, 5096507.198, 5096507.185, 5096507.258, 5096507.293,
              5096507.274, 5096507.247, 5096507.265, 5096507.256, 5096507.28,
              5096507.242, 5096507.212, 5096507.299, 5096507.247, 5096507.273,
              5096507.296, 5096507.293, 5096507.276, 5096507.299, 5096507.285,
              5096507.276, 5096507.217, 5096507.22, 5096507.229, 5096507.307,
              5096507.31, 5096507.328, 5096507.332, 5096507.316, 5096507.315,
              5096507.325, 5096507.309, 5096507.319, 5096507.313, 5096507.304,
              5096507.302, 5096507.058, 5096507.094, 5096507.045, 5096507.072,
              5096507.06, 5096507.107, 5096507.131, 5096507.108, 5096507.105,
              5096507.109, 5096507.101, 5096507.164, 5096507.198, 5096507.159,
              5096507.225, 5096507.212, 5096507.285, 5096507.21, 5096507.219,
              5096507.226, 5096507.214, 5096507.271, 152.553, 152.513, 152.591,
              152.577, 152.561, 152.555, 152.538, 152.548, 152.575, 152.549,
              152.565, 152.55, 152.565, 152.572, 152.553, 152.557, 152.576,
              152.553, 152.557, 152.575, 152.555, 152.579, 152.546, 152.592,
              152.569, 152.574, 152.591, 152.558, 152.571, 152.552, 152.536,
              152.588, 152.555, 152.559, 152.587, 152.559, 152.555, 152.563,
              152.567, 152.558, 152.529, 152.532, 152.548, 152.53, 152.522,
              152.554, 152.541, 152.515, 152.529, 152.519, 152.511, 152.533,
              152.529, 152.55, 152.527, 152.542, 152.534, 152.546, 152.551,
              152.562, 152.538, 152.529, 152.55, 152.535, 152.538, 152.528,
              152.548, 152.558, 152.532, 152.546, 152.522, 152.547, 152.537,
              152.559, 152.567, 152.572, 152.539, 152.542, 152.55, 152.558,
              152.566, 152.535, 152.547, 152.555, 152.545, 152.565, 152.546,
              152.521, 152.549, 152.543, 152.556, 152.527, 152.547, 152.552,
              152.55, 152.543, 152.522, 152.524, 152.532, 152.543, 152.599,
              152.568, 152.587, 152.585, 152.57, 152.556, 152.589, 152.575,
              152.579, 152.584, 152.553, 152.577, 152.58, 152.572, 152.592,
              152.571, 152.531, 152.536, 152.536, 152.551, 152.552, 152.545,
              152.551, 152.521, 152.55, 152.531, 152.545, 152.599, 152.583,
              152.591, 152.589, 152.598, 152.582, 152.59, 152.585, 152.543,
              152.561, 152.549, 152.559, 152.519, 152.554, 152.535, 152.552,
              152.531, 152.552, 152.583, 152.598, 152.576, 152.574, 152.597,
              152.589, 152.586, 152.593, 152.598, 152.59, 152.574, 152.591,
              152.583, 152.543, 152.552, 152.539, 152.579, 152.546, 152.56,
              152.554, 152.525, 152.587, 152.55, 152.585, 152.593, 152.568,
              152.571, 152.58, 152.557, 152.534, 152.541, 152.54, 152.533,
              152.606, 152.603, 152.602, 152.607, 152.609, 152.615, 152.612,
              152.604, 152.601, 152.62, 152.625, 152.625, 152.612, 152.602,
              152.611, 152.606, 152.601, 152.619, 152.601, 152.603, 152.61,
              152.614), dim = c(200L, 3L))

  res <- fit(pt, tolerance = 0.03)

  expect_equal(res$center_x, 312335.577, tolerance = 0.05)
  expect_equal(res$center_y, 5096507.1, tolerance = 0.05)
  expect_equal(res$radius, 0.19, tolerance = 0.03)
  expect_equal(res$covered_arc_degree, 350)
  expect_equal(res$shape_type, "circle")


  if (disp) show(pt, res)
})

test_that("fitting handles complex_slice0", {
  f <- system.file("extdata", "complex_slice0.las", package = "arbor")
  las <- lidR::readLAS(f)
  xyz <- sf::st_coordinates(las)
  res <- fit(xyz, tolerance = 0.03)
  expect_true(is.list(res))


  if (disp) show(xyz, res)
})

test_that("fitting handles complex_slice1", {
  f <- system.file("extdata", "complex_slice1.las", package = "arbor")
  las <- lidR::readLAS(f)
  xyz <- sf::st_coordinates(las)
  res <- fit(xyz, tolerance = 0.03)
  expect_true(is.list(res))

  if (disp) show(xyz, res)

  res <- fit(xyz, tolerance = 0.03, complexity = 3)
  expect_true(is.list(res))

  if (disp) show(xyz, res)
})

test_that("fitting handles complex_slice2", {
  f <- system.file("extdata", "complex_slice2.las", package = "arbor")
  las <- lidR::readLAS(f)
  xyz <- sf::st_coordinates(las)
  res <- fit(xyz, tolerance = 0.03)

  if (disp) show(xyz, res)

  expect_true(is.list(res))
  res <- fit(xyz, tolerance = 0.03, complexity = 3)
  expect_true(is.list(res))

  if (disp) show(xyz, res)
})

test_that("fitting handles complex_slice3", {
  f <- system.file("extdata", "complex_slice3.las", package = "arbor")
  las <- lidR::readLAS(f)
  xyz <- sf::st_coordinates(las)
  res <- fit(xyz, tolerance = 0.03, complexity = 3)
  expect_true(is.list(res))

  if (disp) show(xyz, res)
})

test_that("fitting handles slice_buttress", {
  f <- system.file("extdata", "slice_buttress.las", package = "arbor")
  las <- lidR::readLAS(f)
  xyz <- sf::st_coordinates(las)
  res <- fit(xyz, tolerance = 0.03, complexity = 3)
  expect_true(is.list(res))

  if (disp) show(xyz, res)
})

test_that("fitting handles double_ring_slice0", {
  f <- system.file("extdata", "double_ring_slice0.las", package = "arbor")
  las <- lidR::readLAS(f)
  xyz <- sf::st_coordinates(las)
  res <- fit(xyz, tolerance = 0.03)
  expect_true(is.list(res))

  if (disp) show(xyz, res)

  res <- fit(xyz, tolerance = 0.03, complexity = 3)
  expect_true(is.list(res))

  if (disp) show(xyz, res)
})

test_that("fitting handles double_ring_slice1", {
  f <- system.file("extdata", "double_ring_slice1.las", package = "arbor")
  las <- lidR::readLAS(f)
  xyz <- sf::st_coordinates(las)
  res <- fit(xyz, tolerance = 0.03)
  expect_true(is.list(res))

  if (disp) show(xyz, res)

  res <- fit(xyz, tolerance = 0.03, complexity = 3)
  expect_true(is.list(res))

  if (disp) show(xyz, res)
})

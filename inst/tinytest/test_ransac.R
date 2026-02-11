# Test case 7: RANSAC on noisy circle points
set.seed(123)
n = 500
theta <- seq(0, 2 * pi, length.out = n)
circle_points <- matrix(c(12 + 2 * cos(theta) + rnorm(n, 0, 0.1),
                          18.5 + 2 * sin(theta) + rnorm(n, 0, 0.1),
                          rep(0, n)), ncol = 3, byrow = FALSE)

result7 <- arbor:::ransac_circle(circle_points, num_iterations = 50, inlier_threshold = 0.2)

expect_equal(result7$center_x, 12, tolerance = 0.005)
expect_equal(result7$center_y, 18.5, tolerance = 0.005)
expect_equal(result7$radius, 2, tolerance = 0.02)
expect_equal(result7$covered_arc_degree, 360)


# Test case 8: RANSAC on noisy circle points
set.seed(123)
n = 1000
theta <- seq(0, 2 * pi, length.out = n)
circle_points <- matrix(c(500 + 0.12 * cos(theta) + runif(n, 0, 0.01),
                          500 + 0.12 * sin(theta) + runif(n, 0, 0.01),
                          rep(0, n)), ncol = 3, byrow = FALSE)
rm = circle_points[,1] > 500
circle_points = circle_points[rm,]

result8 <- arbor:::ransac_circle(circle_points, num_iterations = 400)

expect_equal(result8$center_x, 500, tolerance = 0.005)
expect_equal(result8$center_y, 500, tolerance = 0.005)
expect_equal(result8$radius, 0.12, tolerance = 0.05)
expect_equal(result8$covered_arc_degree, 183, tolerance = 3)
expect_equal(result8$percentage_inlier, 1, tolerance = 0.01)



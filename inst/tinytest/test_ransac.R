# Test case 1: Perfect circle
points1 <- matrix(c(0, 1, 0,
                    1, 0, 0,
                    -1, 0, 0),
                  ncol = 3, byrow = TRUE)
result1 <- lidRtls:::fit_circle_on_3_points(points1)
result1 = round(result1, 3)
expect_equal(result1, c(0, 0, 1))

# Test case 2: Points forming an equilateral triangle
points2 <- matrix(c(0, 2, 0,   # Point 1
                    sqrt(3), -1, 0,  # Point 2
                    -sqrt(3), -1, 0), # Point 3
                  ncol = 3, byrow = TRUE)
result2 <- lidRtls:::fit_circle_on_3_points(points2)
result2 = round(result2, 3)
expect_equal(result2, c(0, 0, 2))

# Test case 3: Degenerate case (collinear points)
points3 <- matrix(c(0, 0, 0,   # Point 1
                    1, 1, 0,   # Point 2
                    2, 2, 0),  # Point 3
                  ncol = 3, byrow = TRUE)
result3 <- lidRtls:::fit_circle_on_3_points(points3)
result3 = round(result3, 3)
expect_equal(result3, c(0, 0, 0), info = "Collinear points should return (0,0,0)")

# Test case 4: Larger circle
points4 <- matrix(c(5, 5, 0,   # Point 1
                    10, 5, 0,  # Point 2
                    5, 10, 0), # Point 3
                  ncol = 3, byrow = TRUE)
result4 <- lidRtls:::fit_circle_on_3_points(points4)
result4 = round(result4, 3)
expect_equal(result4, c(7.5,7.5, 3.536))

# Test case 5: Circle with center at (2,2)
points5 <- matrix(c(2, 5, 0,   # Point 1
                    5, 2, 0,   # Point 2
                    -1, 2, 0), # Point 3
                  ncol = 3, byrow = TRUE)
result5 <- lidRtls:::fit_circle_on_3_points(points5)
result5 = round(result5, 3)
expect_equal(result5, c(2,2,3))

# Test case 6: Small circle
points6 <- matrix(c(0, 0.1, 0,   # Point 1
                    0.1, 0, 0,   # Point 2
                    -0.1, 0, 0),  # Point 3
                  ncol = 3, byrow = TRUE)
result6 <- lidRtls:::fit_circle_on_3_points(points6)
result6 = round(result6, 3)
expect_equal(result6, c(0,0,0.1))

# Test case 7: RANSAC on noisy circle points
set.seed(123)
n = 500
theta <- seq(0, 2 * pi, length.out = n)
circle_points <- matrix(c(12 + 2 * cos(theta) + rnorm(n, 0, 0.1),
                          18.5 + 2 * sin(theta) + rnorm(n, 0, 0.1),
                          rep(0, n)), ncol = 3, byrow = FALSE)

result7 <- lidRtls:::ransac_circle(circle_points, num_iterations = 50, inlier_threshold = 0.2)

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

result8 <- lidRtls:::ransac_circle(circle_points, num_iterations = 400)

expect_equal(result8$center_x, 500, tolerance = 0.005)
expect_equal(result8$center_y, 500, tolerance = 0.005)
expect_equal(result8$radius, 0.12, tolerance = 0.05)
expect_equal(result8$covered_arc_degree, 183, tolerance = 3)
expect_equal(result8$percentage_inlier, 1, tolerance = 0.01)


# Test case 8: CPP RANSAC works
set.seed(123)
n = 1000
theta <- seq(0, 2 * pi, length.out = n)
circle_points <- matrix(c(500 + 0.12 * cos(theta) + runif(n, 0, 0.01),
                          500 + 0.12 * sin(theta) + runif(n, 0, 0.01),
                          rep(0, n)), ncol = 3, byrow = FALSE)
rm = circle_points[,1] > 500
circle_points = circle_points[rm,]

result8 <- lidRtls:::ransac_circle_cpp(circle_points, num_iterations = 400, early_exit = 1)

expect_equal(result8$center_x, 500, tolerance = 0.005)
expect_equal(result8$center_y, 500, tolerance = 0.005)
expect_equal(result8$radius, 0.12, tolerance = 0.05)
expect_equal(result8$covered_arc_degree, 183, tolerance = 3)
expect_equal(result8$percentage_inlier, 1, tolerance = 0.01)

# Test case 9: CPP RANSAC works
set.seed(123)
n = 1000
theta <- seq(0, 2 * pi, length.out = n)
circle_points <- matrix(c(500 + 0.12 * cos(theta) + runif(n, 0, 0.01),
                          500 + 0.12 * sin(theta) + runif(n, 0, 0.01),
                          rep(0, n)), ncol = 3, byrow = FALSE)
rm = circle_points[,1] > 500
circle_points = circle_points[rm,]

result9 <- lidRtls:::ransac_circle_cpp(circle_points, num_iterations = 400, early_exit = 0.8)

expect_equal(result9$center_x, 500, tolerance = 0.005)
expect_equal(result9$center_y, 500, tolerance = 0.005)
expect_equal(result9$radius, 0.12, tolerance = 0.05)
expect_equal(result9$covered_arc_degree, 183, tolerance = 3)
expect_true(result9$percentage_inlier > 0.95)



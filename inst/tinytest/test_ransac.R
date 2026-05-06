# Test case 7: RANSAC on noisy circle points
fit = arbor:::fit_circloid_cpp

show = function(pt, res, tol = 0.03)
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

disp = FALSE

set.seed(123)
n = 500
theta <- seq(0, 2 * pi, length.out = n)
circle_points <- matrix(c(12 + 2 * cos(theta) + rnorm(n, 0, 0.1),
                          18.5 + 2 * sin(theta) + rnorm(n, 0, 0.1),
                          rep(0, n)), ncol = 3, byrow = FALSE)

result7 <- fit(circle_points, tolerance = 0.02)

if (disp) show(circle_points, result7)

expect_equal(result7$center_x, 12, tolerance = 0.005)
expect_equal(result7$center_y, 18.5, tolerance = 0.005)
expect_equal(result7$radius, 2, tolerance = 0.025)
expect_equal(result7$covered_arc_degree, 360)


# Test case 8: RANSAC on noisy circle points
set.seed(123)
n = 1000
theta <- seq(0, 2 * pi, length.out = n)
circle_points <- matrix(c(500 + 0.12 * cos(theta) + runif(n, 0, 0.01),
                          500 + 0.12 * sin(theta) + runif(n, 0, 0.01),
                          rep(0, n)), ncol = 3, byrow = FALSE)
rm = circle_points[,1] > 499.95
circle_points = circle_points[rm,]

result8 <- fit(circle_points, tolerance = 0.01)

if (disp) show(circle_points, result8, 0.02)

expect_equal(result8$center_x, 500, tolerance = 0.005)
expect_equal(result8$center_y, 500, tolerance = 0.005)
expect_equal(result8$radius, 0.12, tolerance = 0.05)
expect_equal(result8$covered_arc_degree, 183, tolerance = 3)
expect_equal(result8$percentage_inlier, 100, tolerance = 0.02)


x = structure(
  c(305610.02215, 305610.04015, 305610.05725, 305610.06855,
    305610.05065, 305610.06075, 305609.99775, 305609.98925, 305609.99825,
    305609.99845, 305610.00445, 305610.00355, 305610.04835, 305610.02385,
    305610.01525, 305610.04555, 305610.04175, 305610.03705, 305610.04215,
    305610.02815, 305610.04175, 305610.00565, 305610.01615, 305610.02805,
    305610.03265, 305610.04105, 305610.05135, 305610.03705, 305610.05205,
    305610.05265, 305610.04965, 305610.05565, 305610.06835, 305610.06755,
    305610.04145, 305610.05175, 305610.06255, 305610.05705, 5093250.85515,
    5093250.88565, 5093250.89845, 5093250.89475, 5093250.89655, 5093250.89615,
    5093250.91005, 5093250.91025, 5093250.91315, 5093250.90935, 5093250.93595,
    5093250.92435, 5093250.93635, 5093250.93025, 5093250.91715, 5093250.95825,
    5093250.95265, 5093250.94175, 5093250.94655, 5093250.94845, 5093250.92875,
    5093250.90405, 5093250.94425, 5093250.92905, 5093250.92875, 5093250.90915,
    5093250.90545, 5093250.91385, 5093250.94255, 5093250.94535, 5093250.96135,
    5093250.93885, 5093250.90235, 5093250.91495, 5093250.90015, 5093250.90575,
    5093250.93735, 5093250.92615, 201.7292, 201.7344, 201.7303, 201.7428,
    201.7529, 201.7215, 201.7329, 201.7413, 201.7475, 201.7318, 201.7431,
    201.7212, 201.7089, 201.7371, 201.7557, 201.7268, 201.7549, 201.7448,
    201.7144, 201.7324, 201.7538, 201.7232, 201.7432, 201.7118, 201.7387,
    201.7495, 201.7443, 201.712, 201.7431, 201.7385, 201.7427, 201.706,
    201.7047, 201.7211, 201.7299, 201.71, 201.7077, 201.7411),
  dim = c(38L, 3L), dimnames = list(NULL, c("X", "Y", "Z")))

xc = mean(x[,1])
yc = mean(x[,2])
c = fit(x, tolerance = 0.02)

expect_true(c$center_x - xc < 0.02)
expect_true(c$center_y - yc < 0.02)


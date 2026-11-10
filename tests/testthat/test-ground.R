# A 20 x 20 m plot on a gentle slope: a dense ground, a few trunks and some points far below the ground
set.seed(42)
slope <- function(x) 100 + 0.1 * x

n <- 40000
ground <- data.table::data.table(X = runif(n, 0, 20), Y = runif(n, 0, 20))
ground[, Z := slope(X) + rnorm(n, 0, 0.02)]

trunks <- data.table::rbindlist(lapply(seq(3, 17, by = 7), function(cx) {
  a <- runif(3000, 0, 2 * pi)
  h <- runif(3000, 0, 10)
  data.table::data.table(X = cx + 0.15 * cos(a), Y = 10 + 0.15 * sin(a), Z = slope(cx) + h)
}))

m <- 50
noise <- data.table::data.table(X = runif(m, 1, 19), Y = runif(m, 1, 19))
noise[, Z := slope(X) - runif(m, 1, 3)]

data <- rbind(ground, trunks, noise)
data[, Classification := 0L]
las <- suppressMessages(lidR::LAS(data))
is_noise <- seq_len(nrow(data)) > nrow(ground) + nrow(trunks)

test_that("segment_ground classifies the points below the ground as noise", {
  out <- segment_ground(las)
  cls <- out$Classification + 0L # force a copy: the column was written by C++
  expect_true(all(cls[is_noise] == 7L))
  expect_false(any(cls[!is_noise] == 7L))
  expect_true(sum(cls == 2L) > 0)
  expect_true(all(abs(out$hag[!is_noise]) < 11))
})

test_that("the noise pre-classification can be disabled", {
  data2 <- data.table::copy(data)
  data2[, Classification := 0L] # the previous test wrote the classes in place
  las2 <- suppressMessages(lidR::LAS(data2))
  params <- arbor_parameters_default
  params$ground$noise_enabled <- FALSE
  out <- segment_ground(las2, params)
  cls <- out$Classification + 0L
  expect_false(any(cls == 7L))
})

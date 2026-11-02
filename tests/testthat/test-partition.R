# ScopedPartition (C++) hides the points below the cut height during a stage and puts every point back
# at its original index afterwards. C_scoped_partition() reports what the algorithms see within the
# scope and writes passage = 1 for the kept points, passage = 2 for the hidden ones.

make_cloud <- function(n = 20000L, seed = 42)
{
  set.seed(seed)
  data.table::data.table(
    X = runif(n, 0, 20),
    Y = runif(n, 0, 20),
    Z = runif(n, 0, 10),
    hag = runif(n, -0.1, 2),
    PID = seq_len(n),
    passage = 0L,
    other = rnorm(n))
}

test_that("ScopedPartition hides the low points within the scope only", {
  cloud <- make_cloud()
  ref   <- data.table::copy(cloud)
  cut   <- 0.25

  r <- C_scoped_partition(cloud, cut, FALSE)

  # Within the scope
  expect_equal(r$kept, sum(ref$hag > cut))
  expect_equal(r$total, nrow(ref))
  expect_true(r$kept_above)
  expect_true(r$hidden_below)
  expect_true(r$knn_kept_only)

  # After the scope: size, order and attributes back, KD-tree rebuilt on every point
  expect_equal(r$size_after, nrow(ref))
  expect_true(r$knn_whole)
  expect_identical(cloud[, !"passage"], ref[, !"passage"])

  # The attributes written within the scope followed the points back to their original index
  expect_identical(cloud$passage, ifelse(ref$hag > cut, 1L, 2L))
})

test_that("ScopedPartition restores the point cloud when an error leaves the scope", {
  cloud <- make_cloud()
  ref   <- data.table::copy(cloud)
  cut   <- 0.25

  expect_error(C_scoped_partition(cloud, cut, TRUE), "failure within the scope")

  expect_identical(cloud[, !"passage"], ref[, !"passage"])
  expect_identical(cloud$passage, ifelse(ref$hag > cut, 1L, 2L))
})

test_that("ScopedPartition handles edge cases", {
  # All points kept
  cloud <- make_cloud(1000L)
  ref   <- data.table::copy(cloud)
  r <- C_scoped_partition(cloud, -1, FALSE)
  expect_equal(r$kept, 1000)
  expect_equal(r$size_after, 1000)
  expect_true(r$knn_whole)
  expect_identical(cloud[, !"passage"], ref[, !"passage"])
  expect_all_equal(cloud$passage, 1L)

  # All points hidden
  cloud <- make_cloud(1000L)
  r <- C_scoped_partition(cloud, 10, FALSE)
  expect_equal(r$kept, 0)
  expect_equal(r$size_after, 1000)
  expect_true(r$knn_whole)
  expect_identical(cloud[, !"passage"], ref[, !"passage"])
  expect_all_equal(cloud$passage, 2L)

  # Already partitioned (kept points first, no swap) and reversed (hidden points first, many swaps)
  for (decreasing in c(TRUE, FALSE))
  {
    cloud <- make_cloud(1000L)
    data.table::setorderv(cloud, "hag", order = if (decreasing) -1L else 1L)
    ref <- data.table::copy(cloud)
    r <- C_scoped_partition(cloud, 0.25, FALSE)
    expect_equal(r$kept, sum(ref$hag > 0.25))
    expect_true(r$kept_above && r$hidden_below && r$knn_kept_only && r$knn_whole)
    expect_identical(cloud[, !"passage"], ref[, !"passage"])
    expect_identical(cloud$passage, ifelse(ref$hag > 0.25, 1L, 2L))
  }

  # Empty point cloud
  cloud <- make_cloud(0L)
  r <- C_scoped_partition(cloud, 0.25, FALSE)
  expect_equal(r$kept, 0)
  expect_equal(r$size_after, 0)
})

# The stages that partition the point cloud keep the order and can be run again in place

f   <- system.file("extdata", "9x9.laz", package = "arbor")
las <- lidR::readTLS(f, select = 'xyz', filter = "-keep_random_fraction 0.3")

sink(tempfile())
las <- segment_ground(las)
las <- wood_likelihood(las)
las@data$PID <- seq_len(lidR::npoints(las))
ref <- data.table::copy(las@data)

# The stages write some columns in place (e.g. segment_instance() writes foliage), which R may share between
# LAS objects: each run gets its own copy
fresh <- function(x) { x@data <- data.table::copy(x@data); x }

params <- arbor_parameters_default
params_high <- params
params_high$global$cut_above_ground <- 0.5

first  <- segment_semantic(fresh(las), params)
# A run with another cut height, then in place with the initial one: as a first run
again  <- segment_semantic(segment_semantic(fresh(las), params_high), params)
seeds  <- find_seeds(first)
inst1  <- segment_instance(fresh(first), seeds, params)
# Directly in C++ (the R wrapper resets treeID itself), from tree IDs of a previous run
dirty  <- fresh(first)
dirty@data$treeID <- 99L
segment_instance_cpp(dirty@data, seeds@data, params)
sink()

test_that("segment_semantic() keeps the order and starts from scratch", {
  cut <- params$global$cut_above_ground
  expect_identical(first$PID, ref$PID)
  expect_identical(first@data[, .(X, Y, Z, hag, pwood)], ref[, .(X, Y, Z, hag, pwood)])
  expect_identical(first$UserData == ARBORLOW, !(first$hag > cut))

  expect_identical(again$UserData, first$UserData)
  expect_identical(again$foliage, first$foliage)
  expect_identical(again$passage, first$passage)
})

test_that("segment_instance() keeps the order and starts from scratch", {
  cut <- params$global$cut_above_ground
  expect_identical(inst1$PID, ref$PID)
  expect_identical(inst1@data[, .(X, Y, Z, hag, pwood)], ref[, .(X, Y, Z, hag, pwood)])
  expect_all_equal(inst1$treeID[!(inst1$hag > cut)], 0L)

  expect_identical(dirty$treeID, inst1$treeID)
})

test_that("wood_likelihood() keeps the order and skips the points below the cut height", {
  cut  <- params$global$cut_above_ground
  # Without height above ground, every point is computed
  full <- wood_likelihood(lidR::remove_lasattribute(las, "hag"))$pwood

  expect_identical(las$PID, ref$PID)
  expect_all_equal(las$pwood[!(las$hag > cut)], 0)

  # Far enough above the cut, the neighbourhoods are the same as with every point, except for a few points
  # whose k-th and (k+1)-th neighbours are tied (quantized coordinates): the KD-tree, built on another set
  # of points, may keep the other one, and the difference is platform dependent
  high <- las$hag > cut + 1
  diff <- abs(las$pwood[high] - full[high])
  expect_lt(mean(diff > 1e-6), 1e-4)
  expect_lt(max(diff), 0.1)
})

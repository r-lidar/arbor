f <- system.file("extdata", "9x9.laz", package = "arbor")
las <- lidR::readLAS(f, select = 'xyz', filter = "-keep_random_fraction 0.6")

las <- hybrid_homogeneization(las)
las <- segment_ground(las)
las <- wood_likelihood(las)


test_that("Wood likelihood is valid", {

  expect_true("pwood" %in% names(las))
  expect_true(is.double(las$pwood))
})

las <- segment_semantic(las)

test_that("Semantic segmenation is valid", {
  expect_true("foliage" %in% names(las))
  expect_true(is.integer(las$foliage))

  expect_true("UserData" %in% names(las))
  expect_all_equal(las@data[hag < 0.25]$UserData, ARBORLOW)

  # partitioning. Low points are at the end
  rle = rle(las$UserData)
  expect_equal(rle$values, c(ARBORTREE,ARBORLOW))
})

see <- find_seeds(las)

test_that("Seeds segmenation is valid", {
  expect_is(see, "LAS")
  expect_true("treeID" %in% names(see))
  expect_true(is.integer(see$treeID))
})

las <- segment_instance(las, see)

test_that("Instance segmenation is valid", {
  expect_true("treeID" %in% names(las))
  expect_true(is.integer(las$treeID))

  # partitioning is preserved
  rle = rle(las$UserData)
  expect_equal(rle$values, c(0L,1L))
})

las <- flag_buffer(las, see, -0.75)
las <- flag_small_trees(las, 3)

test_that("Flag buffer is valid", {
  u = sort(unique(las$UserData))
  expect_equal(u, c(ARBORTREE,ARBORLOW,ARBORUNDERSTORY,ARBORBUFFER))
})

qsf <- qsf(las)

test_that("QSF is valid", {
  expect_is(qsf, "qsf")
  expect_is(qsf, "list")
  expect_equal(length(qsf), 22L, tolerance = 1)
})

#plot_semantic(las)
#plot_instance(las)
#lidR::plot(las, color = "UserData")
#plot(qsf, pal = "chocolate4")

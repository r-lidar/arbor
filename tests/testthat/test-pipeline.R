f <- system.file("extdata", "9x9.laz", package = "arbor")
las <- read_tls(f, select = 'xyz', filter = "-keep_random_fraction 0.6")

sink(tempfile())

las <- hybrid_homogeneization(las)
las <- segment_ground(las)
las <- wood_likelihood(las)

# To check that the stages never reorder the points
las@data$PID = 1:lidR::npoints(las)
xyz <- las@data[, .(X, Y, Z, hag)]

test_that("Wood likelihood is valid", {

  expect_true("pwood" %in% names(las))
  expect_true(is.double(las$pwood))
  expect_equal(unique(diff(las$PID)), 1)
})

las <- segment_semantic(las)

#las@data$PID2 = 1:npoints(las)

test_that("Semantic segmenation is valid", {
  expect_true("foliage" %in% names(las))
  expect_true(is.integer(las$foliage))

  expect_true("UserData" %in% names(las))
  expect_all_equal(las@data[hag < 0.25]$UserData, ARBORLOW)
  expect_identical(las$UserData == ARBORLOW, !(las$hag > 0.25))

  # The low points were hidden (partition) during the segmentation only: same order, same coordinates
  expect_identical(las$PID, seq_len(lidR::npoints(las)))
  expect_identical(las@data[, .(X, Y, Z, hag)], xyz)
})

see <- find_seeds(las)

test_that("Seeds segmenation is valid", {
  expect_is(see, "LAS")
  expect_true("treeID" %in% names(see))
  expect_true(is.integer(see$treeID))
})

las <- segment_instance(las, see)

test_that("Overfitting solver works", {

  expect_error(resolve_oversegmentation(las), NA)
})

test_that("Instance segmenation is valid", {
  expect_true("treeID" %in% names(las))
  expect_true(is.integer(las$treeID))

  # treeID is either 0 (no tree) or a valid ID > 0
  expect_false(anyNA(las$treeID))
  expect_true(all(las$treeID >= 0L))

  # Same order, low points untouched by the semantic segmentation flags and without tree
  expect_identical(las$PID, seq_len(lidR::npoints(las)))
  expect_identical(las@data[, .(X, Y, Z, hag)], xyz)
  expect_identical(las$UserData == ARBORLOW, !(las$hag > 0.25))
  expect_all_equal(las$treeID[las$UserData == ARBORLOW], 0L)
})

n   <- lidR::npoints(las)
tid <- las$treeID
low <- las$UserData == ARBORLOW

las <- flag_buffer(las, see, -0.75)
las <- flag_small_trees(las, 1)

test_that("Flag buffer is valid", {
  # Flags only update UserData and never flag LOW points. BUFFER overrides UNDERSTORY
  expect_equal(lidR::npoints(las), n)
  expect_identical(las$treeID, tid)
  expect_all_equal(las$UserData[low], ARBORLOW)
  expect_identical(las$UserData == ARBORBUFFER, flag_buffer(las, see, -0.75)$UserData == ARBORBUFFER)

  u = sort(unique(las$UserData))
  expect_equal(u, c(ARBORTREE,ARBORLOW,ARBORUNDERSTORY,ARBORBUFFER))
})

qsf <- qsf(las)

test_that("QSF is valid", {
  expect_is(qsf, "qsf")
  expect_is(qsf, "list")
  expect_equal(length(qsf), 22L, tolerance = 1)
})

test_that("The segmented point cloud is written and read back with the same types and values", {
  out <- tempfile(fileext = ".laz")
  lidR::writeLAS(las, out)

  # Extra bytes as declared by arbor: no scale nor offset, otherwise rlas reads them as double
  eb <- rlas::read.lasheader(out)[["Variable Length Records"]][["Extra_Bytes"]][["Extra Bytes Description"]]
  types <- c(hag = 10L, pwood = 9L, foliage = 2L, passage = 6L, treeID = 6L) # double, float, char, int, int
  for (name in names(types))
  {
    expect_equal(eb[[name]][["data_type"]], types[[name]], label = paste(name, "data_type"))
    expect_equal(bitwAnd(eb[[name]][["options"]], 24L), 0L, label = paste(name, "scale/offset flags"))
  }

  # Read back with rlas directly: no type conversion by read_tls()
  back <- rlas::read.las(out)
  for (name in c("hag", "pwood", "foliage", "passage", "treeID", "UserData"))
    expect_identical(typeof(back[[name]]), typeof(las@data[[name]]), label = paste(name, "type"))

  expect_equal(nrow(back), lidR::npoints(las))
  for (name in c("foliage", "passage", "treeID", "UserData"))
    expect_identical(back[[name]], las@data[[name]], label = name)
  expect_equal(back$hag, las$hag, tolerance = 1e-6)
  expect_equal(back$pwood, las$pwood, tolerance = 1e-6)
})

sink()

#plot_semantic(las)
#plot_instance(las)
#lidR::plot(las, color = "UserData")
#plot(qsf, pal = "chocolate4")


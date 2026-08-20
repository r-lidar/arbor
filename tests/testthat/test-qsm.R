f <- system.file("extdata", "tree_qsm.laz", package = "arbor")
tree <- lidR::readLAS(f)

sink(tempfile())
qsm1 <- qsm(tree)
sink()

test_that("qsm attributes exist", {
  expect_equal(attr(qsm1, "id"), 0)
  expect_equal(attr(qsm1, "name"), "")
  expect_equal(sf::st_crs(qsm1)$Name, "NAD83 / MTM zone 7")
})

test_that("qsm computes correct dbh", {
  dbh1.3 <- qsm_dbh(qsm1, bh = 1.3)
  dbh2.0 <- qsm_dbh(qsm1, bh = 2.0)

  expect_equal(dbh1.3$dbh, 0.2228, tolerance = 0.005)
  expect_equal(dbh2.0$dbh, 0.2093, tolerance = 0.005)
})

test_that("calling C++ preserves qsm attributes", {
  attr(qsm1, "id") <- 2L
  attr(qsm1, "name") <- "BOUJ45"
  qsm_merchantable(qsm1)
  expect_equal(attr(qsm1, "id"), 2L)
  expect_equal(attr(qsm1, "name"), "BOUJ45")
})

test_that("write and read preserve qsm", {
  attr(qsm1, "id") <- 2L
  attr(qsm1, "name") <- "BOUJ45"

  fqsm <- tempfile(fileext = ".qsm")
  fcsv <- tempfile(fileext = ".csv")
  qsm_write(qsm1, fqsm)
  qsm_write(qsm1, fcsv)
  qsm2 <- qsm_read(fqsm)
  qsm3 <- qsm_read(fcsv)
  v1 <- qsm_volume(qsm1)
  v2 <- qsm_volume(qsm2)
  v3 <- qsm_volume(qsm3)

  expect_equal(v1, v2)
  expect_equal(v1, v3, tolerance = 0.002)
  expect_equal(qsm1$radius, qsm2$radius)
  expect_equal(qsm1$branch_order, qsm2$branch_order)
  expect_equal(sum(qsm1$branch_order == 3), 823L)
  expect_equal(sum(qsm3$branch_order == 3), 823L)
  expect_equal(sum(qsm1$quality == 5), 93L)
  expect_equal(sum(qsm3$quality == 5), 93L)
  expect_equal(attr(qsm2, "id"), 2L)
  expect_equal(attr(qsm2, "name"), "BOUJ45")
  expect_equal(sf::st_crs(qsm2)$Name, "NAD83 / MTM zone 7")
})

test_that("invalid qsm throws error on write", {
  fqsm <- tempfile(fileext = ".qsm")
  hqsm <- qsm1[qsm1$quality > 3, ]
  expect_error(qsm_write(hqsm, fqsm), "Graph is disconnected")
})

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

  expect_equal(dbh1.3$dbh, 0.22, tolerance = 0.005)
  expect_equal(dbh2.0$dbh, 0.226, tolerance = 0.005)
})

test_that("qsm computes correct volume", {
  v1 <- qsm_volume(qsm1)
  expect_equal(v1, 0.476, tolerance = 0.005)
})

test_that("qsm computes correct edges", {
  n1 <- nrow(qsm1)
  expect_equal(n1, 2353L)
})

test_that("qsm computes correct volume per quality", {
  bo1 = qsm1[qsm1$branch_order == 1,]
  bo3 = qsm1[qsm1$branch_order == 3,]
  bq4 = qsm1[qsm1$quality == 4,]
  bq5 = qsm1[qsm1$quality == 5,]
  m1  = qsm_merchantable(qsm1)

  V   = qsm_volume(qsm1)
  vo1 = qsm_volume(bo1)
  vo3 = qsm_volume(bo3)
  vq4 = qsm_volume(bq4)
  vq5 = qsm_volume(bq5)
  vm1 = qsm_volume(m1)

  expect_equal(vo1/V*100, 81.5, tolerance = 0.0005)
  expect_equal(vo3/V*100, 2.9, tolerance = 0.0005)
  expect_equal(vq4/V*100, 3.06, tolerance = 0.0005)
  expect_equal(vq5/V*100, 56.56, tolerance = 0.0005)
  expect_equal(vm1/V*100, 82.73, tolerance = 0.0005)
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
  expect_equal(sum(qsm3$branch_order == 1), sum(qsm1$branch_order == 1))
  expect_equal(sum(qsm3$branch_order == 3), sum(qsm1$branch_order == 3))
  expect_equal(sum(qsm3$quality == 5), sum(qsm1$quality == 5))
  expect_equal(attr(qsm2, "id"), 2L)
  expect_equal(attr(qsm2, "name"), "BOUJ45")
  expect_equal(sf::st_crs(qsm2)$Name, "NAD83 / MTM zone 7")
})

test_that("invalid qsm throws error on write", {
  fqsm <- tempfile(fileext = ".qsm")
  hqsm <- qsm1[qsm1$quality > 3, ]
  expect_error(qsm_write(hqsm, fqsm), "Graph is disconnected")
})

test_that("invalid qsm generate a placeholder", {
  f <- system.file("extdata", "tree_6307.laz", package = "arbor")
  las = lidR::readLAS(f)
  sink(tempfile())
  q = expect_warning(qsm(las))
  sink()
  expect_equal(nrow(q), 1)
  expect_equal(qsm_dbh(q)$dbh, 0)
  q
})


test_that("Topology issues auto repair attempt", {

  sink(tempfile())
  qsm_multi_root <- data.frame(
    startX = c(0, 0, 1),
    startY = c(0, 0, 0),
    startZ = c(0, 1, 1),
    endX   = c(0, 0, 1),
    endY   = c(0, 0, 0),
    endZ   = c(1, 2, 2),
    cyl_ID    = c(2L, 3L, 5L),   # target node
    parent_ID = c(1L, 2L, 4L),   # source node
    axis_ID      = c(1L, 1L, 2L),
    branch_order = c(0L, 0L, 1L),
    radius         = c(0.20, 0.15, 0.05),
    dist_to_root   = c(1, 2, 1),
    subtree_length = c(2, 1, 1),
    quality        = c(4L, 4L, 4L),
    stringsAsFactors = FALSE
  )

  x = arbor:::as_qsm(qsm_multi_root)
  expect_error(print(x), "Graph validation failed")

  qsm_multi_root_fixed = arbor:::qsm_autorepair_cpp(qsm_multi_root)

  x = arbor:::as_qsm(qsm_multi_root_fixed)
  expect_error(print(x), NA)
  sink(NULL)
})

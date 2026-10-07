f <- system.file("extdata", "tree_qsm.laz", package="arbor")
tree <- lidR::readLAS(f)
qsm = qsm(tree)
qsf = arbor:::as_qsf(list(qsm))

test_that("st_crs() works", {
  expect_s3_class(st_crs(qsm), "crs")
  expect_s3_class(st_crs(qsf), "crs")
})

test_that("CRS can be assigned to a qsm", {
  st_crs(qsm) <- 32734
  st_crs(qsf) <- 32734
  expect_s3_class(qsm, "qsm")
  expect_s3_class(qsf, "qsf")
})

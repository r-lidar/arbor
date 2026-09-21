f <- system.file("extdata", "tree_qsm.laz", package = "arbor")
tree <- lidR::readLAS(f)
qsf_obj <- qsf(tree)

sink(tempfile())

test_that("qsf produces stats with required columns", {
  stats <- qsm_stats(qsf_obj)

  expect_true("treeID" %in% names(stats))
  expect_true("V" %in% names(stats))
  expect_true("H" %in% names(stats))
})

test_that("qsf_write(format = 'qsm') and qsf_read() round-trip a forest via a .qsf manifest", {
  dir <- tempfile("qsf_roundtrip_")
  qsf_write(qsf_obj, dir, format = "qsm")

  manifest <- file.path(dir, paste0(basename(dir), ".qsf"))
  expect_true(file.exists(manifest))

  reloaded <- qsf_read(manifest)

  expect_s3_class(reloaded, "qsf")
  expect_equal(length(reloaded), length(qsf_obj))
  expect_equal(sort(names(reloaded)), sort(names(qsf_obj)))
})

test_that("qsf_write via a .qsf manifest", {
  path <- tempfile(fileext = ".qsf")
  qsf_write(qsf_obj, path)

  expect_true(file.exists(path))

  reloaded <- qsf_read(path)

  expect_s3_class(reloaded, "qsf")
  expect_equal(length(reloaded), length(qsf_obj))
  expect_equal(sort(names(reloaded)), sort(names(qsf_obj)))
})



sink()

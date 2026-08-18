f <- system.file("extdata", "tree_qsm.laz", package = "arbor")
tree <- lidR::readLAS(f)

test_that("qsf produces stats with required columns", {
  qsf_obj <- qsf(tree)
  stats <- qsm_stats(qsf_obj)

  expect_true("treeID" %in% names(stats))
  expect_true("V" %in% names(stats))
  expect_true("H" %in% names(stats))
})

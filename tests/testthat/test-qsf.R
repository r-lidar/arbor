sink(tempfile())

f <- system.file("extdata", "tree_qsm.laz", package = "arbor")
tree <- lidR::readLAS(f)
qsf_obj <- qsf(tree)

sink()

f <- system.file("extdata", "oak-plantation.qsf", package="arbor")
test_qsf <- qsf_read(f)

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

test_that("qsf_log parses coded messages", {
  log <- qsf_log(test_qsf)
  expect_s3_class(log, "qsf_log")
  expect_true(setequal(names(log), c("W2", "W3")))

  l <- log[["W2"]]
  expect_equal(l$code, "W2")
  expect_equal(l$type, "warning")
  expect_equal(l$label, "No valid measure")
  expect_equal(l$treeID, c(149L, 177L, 290L, 487L, 516L, 750L, 852L))

  l <- log[["W3"]]
  expect_equal(l$type, "warning")
  expect_equal(l$treeID, c(102L, 165L, 172L, 175L, 215L, 231L, 232L, 259L, 524L, 659L, 891L))
})

test_that("qsf_log returns an empty list when there are no messages", {
  test_qsf <- arbor:::as_qsf(qsf_obj[1])
  attr(test_qsf[[1]], "message") <- NULL
  expect_equal(qsf_log(test_qsf), list())
})

test_that("qsf_select excludes/keeps trees by log code and type", {
  attr(test_qsf[[3]], "message") <- "[E0] [Error] Boom"
  
  q <- qsf_select(test_qsf, code = "W2")
  expect_equal(length(q), 7L)

  q <- qsf_select(test_qsf, code = "W3", invert = TRUE)
  expect_equal(length(q), 47L)

  q <- qsf_select(test_qsf, type = "error")
  expect_equal(length(q), 1L)

  expect_s3_class(q, "qsf")
})

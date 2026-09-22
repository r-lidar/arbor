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

test_that("qsf_log parses coded messages, legacy uncoded messages and a mix of both", {
  qsm_legacy <- qsf_obj[[1]]
  qsm_coded  <- qsf_obj[[1]]
  qsm_multi  <- qsf_obj[[1]]

  attr(qsm_legacy, "message") <- "[No wood point] This tree has no point labelled as wood"
  attr(qsm_coded, "message")  <- "[W0] [No wood point] This tree has no point labelled as wood"
  attr(qsm_multi, "message")  <- c("[W2] [No valid measure] retry", "[E1] [Fatal] boom")

  test_qsf <- arbor:::as_qsf(list(`10` = qsm_legacy, `20` = qsm_coded, `30` = qsm_multi))

  log <- qsf_log(test_qsf)

  expect_s3_class(log, "qsf_log")
  expect_true(setequal(names(log), c("No wood point", "W0", "W2", "E1")))

  legacy <- log[["No wood point"]]
  expect_true(is.na(legacy$code))
  expect_true(is.na(legacy$type))
  expect_equal(legacy$label, "No wood point")
  expect_equal(legacy$treeID, 10L)

  coded <- log[["W0"]]
  expect_equal(coded$code, "W0")
  expect_equal(coded$type, "warning")
  expect_equal(coded$label, "No wood point")
  expect_equal(coded$treeID, 20L)

  retry <- log[["W2"]]
  expect_equal(retry$type, "warning")
  expect_equal(retry$treeID, 30L)

  err <- log[["E1"]]
  expect_equal(err$code, "E1")
  expect_equal(err$type, "error")
  expect_equal(err$treeID, 30L)
})

test_that("qsf_log returns an empty list when there are no messages", {
  test_qsf <- arbor:::as_qsf(qsf_obj[1])
  attr(test_qsf[[1]], "message") <- NULL
  expect_equal(qsf_log(test_qsf), list())
})

test_that("qsf_filter excludes/keeps trees by log code and type", {
  qsm_ok    <- qsf_obj[[1]]
  qsm_warn  <- qsf_obj[[1]]
  qsm_error <- qsf_obj[[1]]

  attr(qsm_warn, "message")  <- "[W0] [No wood point] This tree has no point labelled as wood"
  attr(qsm_error, "message") <- "[E1] [Fatal] boom"

  test_qsf <- arbor:::as_qsf(list(`1` = qsm_ok, `2` = qsm_warn, `3` = qsm_error))

  excluded <- qsf_filter(test_qsf, code = "W0")
  expect_equal(sort(names(excluded)), c("1", "3"))

  kept <- qsf_filter(test_qsf, code = "W0", invert = TRUE)
  expect_equal(names(kept), "2")

  errors_only <- qsf_filter(test_qsf, type = "error", invert = TRUE)
  expect_equal(names(errors_only), "3")

  expect_s3_class(excluded, "qsf")
})

test_that("qsf_filter filters by dbh and height ranges", {
  dbh_val <- qsm_dbh(qsf_obj[[1]])$dbh
  h_val   <- qsm_height(qsf_obj[[1]])

  test_qsf <- arbor:::as_qsf(qsf_obj[1])

  in_range  <- qsf_filter(test_qsf, dbh = c(dbh_val - 0.01, dbh_val + 0.01))
  out_range <- qsf_filter(test_qsf, dbh = c(dbh_val + 1, dbh_val + 2))

  expect_equal(length(in_range), 1)
  expect_equal(length(out_range), 0)

  in_h  <- qsf_filter(test_qsf, max_height = c(h_val - 1, h_val + 1))
  out_h <- qsf_filter(test_qsf, max_height = c(h_val + 10, h_val + 20))

  expect_equal(length(in_h), 1)
  expect_equal(length(out_h), 0)
})

test_that("qsm_message stays backward compatible with legacy and coded messages", {
  qsm_legacy <- qsf_obj[[1]]
  attr(qsm_legacy, "message") <- "[No wood point] This tree has no point labelled as wood"
  expect_equal(qsm_message(qsm_legacy, short = TRUE), "No wood point")

  qsm_coded <- qsf_obj[[1]]
  attr(qsm_coded, "message") <- "[W0] [No wood point] This tree has no point labelled as wood"
  expect_equal(qsm_message(qsm_coded, short = TRUE), "No wood point")
})



sink()

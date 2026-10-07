f <- system.file("extdata", "oak-plantation.qsf", package="arbor")
test_qsf <- qsf_read(f)

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
  void = qsf_filter_errors(test_qsf)
  expect_length(void, 0)
  expect_s3_class(void, "qsf")
})

test_that("qsf_filter excludes/keeps trees by log code and type", {
  attr(test_qsf[[3]], "message") <- "[E0] [Error] Boom"
  
  q <- qsf_filter(test_qsf, code = "W2")
  expect_equal(length(q), 7L)

  q <- qsf_filter(test_qsf, code = "W3", invert = TRUE)
  expect_equal(length(q), 47L)

  q <- qsf_filter(test_qsf, type = "error")
  expect_equal(length(q), 1L)

  expect_s3_class(q, "qsf")
})

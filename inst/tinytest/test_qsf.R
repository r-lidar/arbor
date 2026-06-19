f <- system.file("extdata", "tree_qsm.laz", package="arbor")
tree <- lidR::readLAS(f)

qsf = qsf(tree)
qsf

stats = qsm_stats(qsf)

expect_true("treeID" %in% names(stats))
expect_true("V" %in% names(stats))
expect_true("H" %in% names(stats))

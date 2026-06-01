f <- system.file("extdata", "tree_qsm.laz", package="arbor")
tree <- lidR::readLAS(f)
qsm1 = qsm(tree)

dbh1.3 = qsm_dbh(qsm1, bh = 1.3)
dbh2.0 = qsm_dbh(qsm1, bh = 2.0)

expect_equal(dbh1.3$dbh,  0.224, tolerance = 0.001)
expect_equal(dbh2.0$dbh,  0.219, tolerance = 0.001)

fqsm = tempfile(fileext = ".qsm")

qsm_write(qsm1, fqsm)
qsm2 = qsm_read(fqsm)

v1 = qsm_volume(qsm1)
v2 = qsm_volume(qsm2)

expect_equal(v1, v2)
expect_equal(qsm1$radius, qsm2$radius)
expect_equal(qsm1$branch_order, qsm2$branch_order)

hqsm = qsm1[qsm1$quality >3,]

expect_error(qsm_write(hqsm, fqsm), "Graph is disconnected")

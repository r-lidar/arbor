f <- system.file("extdata", "tree_qsm.laz", package="arbor")
tree <- lidR::readLAS(f)
qsm1 <- qsm(tree)

# attr exist
expect_equal(attr(qsm1, "id"), 0)
expect_equal(attr(qsm1, "name"), "")


# correct dbh
dbh1.3 = qsm_dbh(qsm1, bh = 1.3)
dbh2.0 = qsm_dbh(qsm1, bh = 2.0)

expect_equal(dbh1.3$dbh,  0.2260, tolerance = 0.005)
expect_equal(dbh2.0$dbh,  0.2200, tolerance = 0.005)

# calling C++ preserve attr
attr(qsm1, "id") = 2L
attr(qsm1, "name") = "BOUJ45"
merch = qsm_merchantable(qsm1)
expect_equal(attr(qsm1, "id"), 2L)
expect_equal(attr(qsm1, "name"), "BOUJ45")


# write preserves qsm
fqsm = tempfile(fileext = ".qsm")
fcsv = tempfile(fileext = ".csv")
qsm_write(qsm1, fqsm)
qsm_write(qsm1, fcsv)
qsm2 = qsm_read(fqsm)
qsm3 = qsm_read(fcsv)
v1 = qsm_volume(qsm1)
v2 = qsm_volume(qsm2)
v3 = qsm_volume(qsm3)

expect_equal(v1, v2)
expect_equal(v1, v3, tolerance = 0.001)
expect_equal(qsm1$radius, qsm2$radius)
expect_equal(qsm1$branch_order, qsm2$branch_order)
expect_equal(sum(qsm1$branch_order == 3), 904L)
expect_equal(sum(qsm3$branch_order == 3), 904L)
expect_equal(sum(qsm1$quality == 3), 176L)
expect_equal(sum(qsm3$quality == 3), 176L)
expect_equal(attr(qsm2, "id"), 2L)
expect_equal(attr(qsm2, "name"), "BOUJ45")

# invalid qsm throw error
hqsm = qsm1[qsm1$quality >3,]
expect_error(qsm_write(hqsm, fqsm), "Graph is disconnected")


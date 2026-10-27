# Inventory merging (port of arborInventoryMerger). The bindings take the point cloud as
# a data.frame and the inventory as a data.frame with ID_Arbre, X, Y, Z, DHP, POM.

quiet <- function(expr) { out <- NULL ; utils::capture.output(out <- expr) ; out }

# ==============================================================================
# Synthetic scene: vertical cylinders standing on a ground at elevation Z0
#
#   pc ID   inventory ID   position     situation
#   1       101            (0, 0)       simple match
#   101     102            (3, 0)       pc ID collides with inventory 101 -> disentangled
#   3       103 and 104    (6, 0)       two stems segmented as a single tree -> conflict
#                          (6.8, 0)     resolved by re-segmentation
#   -       105            (10, 0)      no point -> undetected
#   NA      -              scattered    unassigned points, never remapped
# ==============================================================================

X0 <- 772300 ; Y0 <- 379800 ; Z0 <- 660 # Realistic UTM coordinates

stem <- function(x, y, r, id, n = 4000, h = 5)
{
  z <- stats::runif(n, 0, h)
  a <- stats::runif(n, 0, 2*pi)
  data.frame(
    X = X0 + x + r * cos(a) + stats::rnorm(n, 0, 0.003),
    Y = Y0 + y + r * sin(a) + stats::rnorm(n, 0, 0.003),
    Z = Z0 + z,
    hag = z,
    treeID = as.integer(id),
    foliage = 0L)
}

set.seed(42)
ground <- data.frame(X = X0 + stats::runif(500, -2, 12), Y = Y0 + stats::runif(500, -3, 3), Z = Z0, hag = 0, treeID = NA_integer_, foliage = 1L)
scene <- rbind(stem(0, 0, 0.20, 1), stem(3, 0, 0.15, 101), stem(6, 0, 0.15, 3), stem(6.8, 0, 0.15, 3), ground)
#x = plot_instance(lidR::LAS(scene))

inventory <- data.frame(
  ID_Arbre = c(101L, 102L, 103L, 104L, 105L),
  X   = X0 + c(0, 3, 6, 6.8, 10),
  Y   = Y0 + c(0, 0, 0, 0, 0),
  Z   = Z0 + 1.3,
  DHP = c(0.40, 0.30, 0.30, 0.30, 0.30),
  POM = 1.3)

is_stem <- function(x) abs(scene$X - (X0 + x)) < 0.3 & !is.na(scene$treeID)

test_that("the inventory must be valid", {
  bad <- inventory ; bad$ID_Arbre[1] <- 0L
  expect_error(arbor:::aim_matching_table_cpp(scene, bad), "invalid tree id")

  bad <- inventory ; bad$ID_Arbre[1] <- NA
  expect_error(arbor:::aim_matching_table_cpp(scene, bad), "integral")

  bad <- inventory ; bad$ID_Arbre[1] <- 1.5
  expect_error(arbor:::aim_matching_table_cpp(scene, bad), "integral")

  bad <- inventory ; bad$ID_Arbre[2] <- 101L
  expect_error(arbor:::aim_matching_table_cpp(scene, bad), "already exists")

  bad <- inventory ; bad$DHP[1] <- 0
  expect_error(arbor:::aim_matching_table_cpp(scene, bad), "radius")

  bad <- inventory ; bad$POM[1] <- NA
  expect_error(arbor:::aim_matching_table_cpp(scene, bad), "pom")

  bad <- inventory ; bad$X[1] <- NA
  expect_error(arbor:::aim_matching_table_cpp(scene, bad), "position")

  expect_error(arbor:::aim_matching_table_cpp(scene, inventory[, -6]), "POM")
})

test_that("aim_disentangle_ids_cpp renumbers point cloud IDs that exist in the inventory", {
  original <- scene$treeID
  id <- quiet(arbor:::aim_disentangle_ids_cpp(scene, inventory))

  expect_identical(scene$treeID, original) # input untouched
  expect_true(is.integer(id))
  expect_equal(length(id), nrow(scene))

  # 101 collides with the inventory and becomes max(pc IDs, inventory IDs) + 1 = 106
  expect_true(all(id[original %in% 101L] == 106L))
  expect_false(any(id %in% inventory$ID_Arbre))

  # Other IDs and unassigned points are untouched
  keep <- !original %in% 101L
  expect_identical(id[keep], original[keep])

  # Nothing to do
  expect_identical(quiet(arbor:::aim_disentangle_ids_cpp(scene[!original %in% 101L, ], inventory)), original[keep])
})

scene$treeID <- quiet(arbor:::aim_disentangle_ids_cpp(scene, inventory))

test_that("aim_matching_table_cpp refuses entangled IDs", {
  entangled <- scene
  entangled$treeID[entangled$treeID %in% 106L] <- 101L
  expect_error(quiet(arbor:::aim_matching_table_cpp(entangled, inventory)), "disentangle")
})

test_that("aim_matching_table_cpp matches point cloud trees and inventory trees", {
  mt <- quiet(arbor:::aim_matching_table_cpp(scene, inventory))

  expect_is(mt, "data.frame")
  expect_named(mt, c("treeID", "targetID", "nfit"))
  expect_true(is.integer(mt$treeID))

  mt <- mt[order(mt$targetID, mt$treeID), ]
  expect_equal(mt$targetID, c(101L, 102L, 103L, 104L, 105L))
  expect_equal(mt$treeID,   c(1L, 106L, 3L, 3L, NA))

  # Detected trees are fitted on several slices, the undetected one has 0
  expect_true(all(mt$nfit[1:4] > 1))
  expect_equal(mt$nfit[5], 0L)
})

mt <- quiet(arbor:::aim_matching_table_cpp(scene, inventory))

test_that("aim_fit_tree_cpp returns the fits of a tree", {
  res <- quiet(arbor:::aim_fit_tree_cpp(scene, inventory[1, ]))

  expect_is(res, "list")
  expect_true(length(res) > 1)

  fit <- res[[1]]
  expect_true(fit$success)
  expect_equal(fit$radius, 0.2, tolerance = 0.01)
  expect_equal(fit$center_x, X0, tolerance = 0.01)
  expect_equal(fit$center_y, Y0, tolerance = 0.01)
  expect_true(fit$covered_arc_degree > 180)
  expect_equal(fit$to_be_merge_id, 1L)
  expect_is(fit$nodes, "matrix")
  expect_named(fit$seeds, c("X", "Y", "Z"))
  expect_equal(nrow(fit$seeds), length(fit$inliers))

  # Undetected tree: a single pseudo-fit at the inventory position
  res <- quiet(arbor:::aim_fit_tree_cpp(scene, inventory[5, ]))
  expect_equal(length(res), 1L)
  expect_equal(res[[1]]$center_x, X0 + 10)
  expect_identical(res[[1]]$to_be_merge_id, NA_integer_)
  expect_equal(nrow(res[[1]]$seeds), 1L)

  expect_error(arbor:::aim_fit_tree_cpp(scene, inventory), "multiple trees")
})

test_that("aim_reassign_ids_cpp applies the matching table and leaves conflicts", {
  original <- scene$treeID
  id <- quiet(arbor:::aim_reassign_ids_cpp(scene, mt))

  expect_identical(scene$treeID, original) # input untouched
  expect_true(all(id[original %in% 1L] == 101L))
  expect_true(all(id[original %in% 106L] == 102L))
  expect_true(all(id[original %in% 3L] == 3L))  # conflict 3 -> 103,104 left untouched
  expect_true(all(is.na(id[is.na(original)])))
})

test_that("aim_resolve_multi_matching_cpp splits trees matched to several inventory trees", {
  s <- scene
  s$treeID <- quiet(arbor:::aim_reassign_ids_cpp(s, mt))
  before <- s$treeID
  id <- quiet(arbor:::aim_resolve_multi_matching_cpp(s, mt, inventory))

  expect_identical(s$treeID, before) # input untouched

  # Each stem goes to its own inventory tree. Points below cut_above_ground keep their ID.
  high <- scene$hag > 0.5
  expect_gt(mean(id[is_stem(6) & high] == 103L), 0.99)
  expect_gt(mean(id[is_stem(6.8) & high] == 104L), 0.99)

  # Other trees are untouched
  other <- !(is_stem(6) | is_stem(6.8))
  expect_identical(id[other], s$treeID[other])
})

test_that("aim_merge_inventory_cpp runs the whole pipeline", {
  original <- scene$treeID
  s <- scene
  s$treeID <- quiet(arbor:::aim_reassign_ids_cpp(s, mt))
  expected <- quiet(arbor:::aim_resolve_multi_matching_cpp(s, mt, inventory))

  id <- quiet(arbor:::aim_merge_inventory_cpp(scene, inventory))

  expect_identical(scene$treeID, original) # input untouched
  expect_identical(id, expected)
  expect_true(all(c(101L, 102L, 103L, 104L) %in% id))
  expect_false(105L %in% id)
})

# ==============================================================================
# Reassignment rules on a tiny point cloud (one point per tree)
# ==============================================================================

ids <- function(...) { v <- c(...) ; data.frame(X = as.numeric(seq_along(v)), Y = 0, Z = 0, treeID = as.integer(v)) }
table <- function(tree, target) data.frame(treeID = as.integer(tree), targetID = as.integer(target))

test_that("aim_reassign_ids_cpp never merges trees by accident", {
  # Pre-existing ID equal to a target is moved to a fresh ID (max + n_src + 1 = 3 + 1 + 1)
  expect_identical(quiet(arbor:::aim_reassign_ids_cpp(ids(1, 2, 3), table(1, 2))), c(2L, 5L, 3L))

  # Swap
  expect_identical(quiet(arbor:::aim_reassign_ids_cpp(ids(1, 2, 3), table(c(1, 2), c(2, 1)))), c(2L, 1L, 3L))

  # Identity rows are ignored and do not hide a collision: 3 -> 5 must not merge with 5
  expect_identical(quiet(arbor:::aim_reassign_ids_cpp(ids(3, 5), table(c(5, 3), c(5, 5)))), c(5L, 7L))

  # Duplicated treeID: not reassigned. Other rules are applied.
  expect_identical(quiet(arbor:::aim_reassign_ids_cpp(ids(1, 2), table(c(1, 1, 2), c(10, 11, 12)))), c(1L, 12L))

  # NA and 0 are never remapped, NA rows of the table are ignored
  expect_identical(quiet(arbor:::aim_reassign_ids_cpp(ids(NA, 0, 1), table(c(NA, 1), c(9, 10)))), c(NA, 0L, 10L))

  expect_error(arbor:::aim_reassign_ids_cpp(ids(1), table(1, NA)), "NA")
})

# ==============================================================================
# Local only: regression on a real plot (Mathilde quadrat 2, 36 M points). The
# data are too large to be shipped with the package.
# ==============================================================================

test_that("Inventory merging on a real plot",
{
  dir  <- "/home/jr/Documents/r-lidar inc/clients/ird/Data_4_JRR/Quart_ha_&_inventaire/Mathilde"

  skip_on_cran()
  skip_on_ci()
  skip_if_not(dir.exists(dir))

  flas <- file.path(dir, "quadrat_2_buf15m_decimated25_output/quadrat_2_buf15m_decimated25_segmented.laz")
  fdtm <- file.path(dir, "quadrat_2_buf15m_decimated25_output/quadrat_2_buf15m_decimated25_dtm.tif")
  finv <- file.path(dir, "inventaire.txt")
  skip_if_not(all(file.exists(c(flas, fdtm, finv))), "Local test data not available")
  skip_if_not_installed("terra")

  las <- lidR::readTLS(flas, select = "xyz0")
  data <- las@data[, c("X", "Y", "Z", "hag", "treeID", "foliage")]
  rm(las)

  dtm <- terra::rast(fdtm)
  raw <- utils::read.csv(finv)
  gnd <- terra::extract(dtm, cbind(raw$CenterX, raw$CenterY))[[1]]
  inv <- data.frame(ID_Arbre = raw$IDtree, X = raw$CenterX, Y = raw$CenterY, Z = raw$Z, DHP = raw$DBH, POM = raw$Z - gnd)

  data$treeID <- quiet(arbor:::aim_disentangle_ids_cpp(data, inv))
  expect_false(any(data$treeID %in% inv$ID_Arbre))

  mt <- quiet(arbor:::aim_matching_table_cpp(data, inv))
  mt <- mt[order(mt$targetID, mt$treeID), ]

  expect_equal(mt$targetID, c(21L, 24L, 33L, 34L, 35L, 36L, 36L, 37L, 38L, 39L, 40L, 40L, 41L, 42L, 43L, 44L, 45L, 45L, 47L, 48L, 49L, 50L, 54L, 55L, 56L, 57L, 58L, 59L, 60L, 61L, 199L))
  expect_equal(mt$treeID,   c(28L, 169L, 70L, 70L, 71L, 123L, 893L, 122L, 123L, 123L, 123L, 893L, 5401L, 521L, 92L, 30L, 32L, 582L, 27L, 5L, 16L, 5403L, NA, 385L, 5399L, 5399L, 5399L, 5399L, 5399L, 5398L, 3L))

  id <- quiet(arbor:::aim_merge_inventory_cpp(data, inv))

  s <- data
  s$treeID <- quiet(arbor:::aim_reassign_ids_cpp(s, mt))
  expect_identical(id, quiet(arbor:::aim_resolve_multi_matching_cpp(s, mt, inv)))

  # All inventory trees are in the segmentation except the undetected one
  expect_setequal(intersect(inv$ID_Arbre, id), setdiff(inv$ID_Arbre, 54L))
})

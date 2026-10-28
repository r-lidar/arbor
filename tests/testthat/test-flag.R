# A 30 x 30 m plot on a gentle slope: a 0.5 m grid of LOW ground points and a 7 x 7 grid of trees
# 4.5 m apart whose expected flags are known. Trees are a trunk of radius 15 cm with a crown, made
# of rings of k points every 20 cm. To see a result in 3D: lidR::plot(out, color = "UserData")

slope <- function(x) 0.02 * x

make_tree <- function(id, x, y, h0 = 0, h1 = 15, k = 20L, crown = 1.2, lean = 0)
{
  h <- seq(h0, h1, by = 0.2)
  t <- (h - h0) / max(h1 - h0, 1e-6)
  r <- pmax(0.15, crown * sin(pi * pmax(t - 0.4, 0) / 0.6))
  a <- 2 * pi * seq_len(k) / k
  data.table::data.table(
    X = x + as.vector(outer(cos(a), r)) + rep(t, each = k) * lean,
    Y = y + as.vector(outer(sin(a), r)) + rep(t, each = k) * lean,
    Z = slope(x) + rep(h, each = k),
    hag = rep(h, each = k),
    treeID = as.integer(id),
    UserData = ARBORTREE)
}

pos   <- c(1.5, 6, 10.5, 15, 19.5, 24, 28.5)
trees <- data.table::CJ(x = pos, y = pos)
trees[, `:=`(id = .I, h0 = 0, h1 = 8 + (.I %% 5) * 3, k = 20L, lean = 0)]
trees[, edge := pmin(x, 30 - x, y, 30 - y)]                  # distance from the base to the plot edge
trees[, small := x %in% c(6, 24) & y == 15 | x == 15 & y == 6 | x == 1.5 & y == 6] # 1.5 m tall, one in the buffer
trees[, fragment := x == 10.5 & y == 19.5]                    # 3 points floating at 3 m: flagged as small (no point below the threshold), no root without seeds
trees[, leaning := x == 24 & y == 24]                         # base 6 m from the edge, top 1.5 m
trees[small == TRUE, h1 := 1.5]
trees[fragment == TRUE, `:=`(h0 = 3, h1 = 3.4, k = 1L)]
trees[leaning == TRUE, lean := 4.5]

g <- seq(0, 30, by = 0.5)
ground <- data.table::CJ(X = g, Y = g)
ground[, `:=`(Z = slope(X), hag = 0, treeID = 0L, UserData = ARBORLOW)]

# Unsegmented shrubs (treeID 0, NA, -1): one in the middle of the plot, one in the buffer
shrub <- rbind(make_tree(0, 12.75, 12.75, h1 = 1, k = 6L, crown = 0.5), make_tree(0, 1, 12.75, h1 = 1, k = 6L, crown = 0.5))
shrub[, treeID := rep(c(0L, NA_integer_, -1L), length.out = .N)]

scene <- data.table::rbindlist(c(
  list(ground, shrub),
  lapply(seq_len(nrow(trees)), function(i) with(trees[i], make_tree(id, x, y, h0, h1, k, lean = lean)))))
scene <- suppressMessages(lidR::LAS(scene, check = FALSE))

is_ground <- scene$UserData == ARBORLOW
is_shrub  <- !is_ground & !(scene$treeID %in% trees$id)
in_shrub  <- is_shrub & scene$X > 5 # the shrub in the middle

# Seeds at the base of each tree
seeds <- suppressMessages(lidR::LAS(trees[, list(X = x, Y = y, Z = 0, treeID = id)], check = FALSE))

# Check that each tree has a single UserData: 'value' if 'flagged', ARBORTREE otherwise
expect_flags <- function(out, flagged, value)
{
  d <- out@data[treeID %in% trees$id, list(ud = unique(UserData)), by = treeID]
  expect_equal(anyDuplicated(d$treeID), 0L)
  expect_equal(d$ud[match(trees$id, d$treeID)], ifelse(flagged, value, ARBORTREE))
}

test_that("flag_small_trees() flags small and floating trees (no point below the threshold)", {
  ud0 <- data.table::copy(scene$UserData)

  out <- flag_small_trees(scene)
  #lidR::plot(out, color = "UserData")
  expect_flags(out, trees$small | trees$fragment, ARBORUNDERSTORY)
  expect_equal(unique(out$UserData[is_ground]), ARBORLOW)
  expect_equal(unique(out$UserData[is_shrub]), ARBORUNDERSTORY)
  expect_identical(out$treeID, scene$treeID)
  expect_identical(scene$UserData, ud0) # input not modified

  # Previous flags are reverted
  out <- flag_small_trees(out, 1)
  #lidR::plot(out, color = "UserData")
  expect_flags(out, trees$fragment, ARBORUNDERSTORY)
})

test_that("flag_buffer() flags trees whose root is close to the edge", {
  ud0 <- data.table::copy(scene$UserData)

  # The root of the leaning tree is the mean of its 100 lowest points, not of all its points.
  # Unsegmented points are flagged according to their own position
  out <- flag_buffer(scene)
  #lidR::plot(out, color = "UserData")
  expect_flags(out, trees$edge < 5, ARBORBUFFER)
  expect_equal(unique(out$UserData[is_ground]), ARBORLOW)
  expect_equal(unique(out$UserData[in_shrub]), ARBORTREE)
  expect_equal(unique(out$UserData[is_shrub & !in_shrub]), ARBORBUFFER)
  expect_identical(out$treeID, scene$treeID)
  expect_identical(scene$UserData, ud0) # input not modified

  # BUFFER overrides UNDERSTORY whatever the order of the flags
  buf <- out$UserData == ARBORBUFFER
  expect_identical(flag_small_trees(out)$UserData == ARBORBUFFER, buf)
  expect_identical(flag_buffer(flag_small_trees(scene))$UserData == ARBORBUFFER, buf)

  # Seeds are the roots. Roots exactly on the edge are inside. Previous flags are reverted
  out <- flag_buffer(scene, buffer = -100)
  out <- flag_buffer(out, seeds, -6)
  #lidR::plot(out, color = "UserData")
  expect_flags(out, trees$edge < 6, ARBORBUFFER)
})

test_that("flag_buffer() with a polygon", {
  # A square that contains the trees of the lower left quarter
  square <- sf::st_polygon(list(rbind(c(0, 0), c(16, 0), c(16, 16), c(0, 16), c(0, 0))))
  square <- sf::st_sfc(square, crs = sf::st_crs(scene))
  out <- flag_buffer(scene, buffer = square)
  #lidR::plot(out, color = "UserData")
  expect_flags(out, trees$x > 16 | trees$y > 16, ARBORBUFFER)

  # A U shape whose notch contains the trees of the central column in the upper half
  u <- sf::st_polygon(list(rbind(c(0, 0), c(30, 0), c(30, 30), c(17, 30), c(17, 12), c(13, 12), c(13, 30), c(0, 30), c(0, 0))))
  u <- sf::st_sf(id = 1, geometry = sf::st_sfc(u, crs = sf::st_crs(scene)))
  out <- flag_buffer(scene, seeds, u)
  #lidR::plot(out, color = "UserData")
  expect_flags(out, trees$x == 15 & trees$y > 12, ARBORBUFFER)
})

test_that("flag_buffer() shrinks the convex hull like sf::st_buffer()", {
  # Many single point trees in an irregular hull with large coordinates
  set.seed(3)
  n <- 10000L
  ang <- runif(n, 0, 2*pi)
  rad <- 30 * sqrt(runif(n)) * (1 + 0.3*cos(3*ang))
  d <- data.table::data.table(X = 500000 + rad*cos(ang), Y = 5e6 + rad*sin(ang), Z = 0, treeID = seq_len(n), UserData = 0L)
  las <- suppressMessages(lidR::LAS(d, check = FALSE))

  pts    <- sf::st_as_sf(d, coords = c("X", "Y"))
  hull   <- sf::st_convex_hull(sf::st_union(pts))
  inside <- lengths(sf::st_intersects(pts, sf::st_buffer(hull, -3.3))) > 0
  out <- flag_buffer(las, las, -3.3)
  expect_equal(out$UserData, ifelse(inside, ARBORTREE, ARBORBUFFER))
})

test_that("flags validate their inputs", {
  sq    <- rbind(c(0, 0), c(5, 0), c(5, 5), c(0, 5), c(0, 0))
  crs   <- sf::st_crs(scene)
  multi <- sf::st_sfc(sf::st_multipolygon(list(list(sq))), crs = crs)
  hole  <- sf::st_sfc(sf::st_polygon(list(sq * 3, sq + 5)), crs = crs)
  other <- sf::st_sfc(sf::st_polygon(list(sq)), crs = 2949)

  expect_error(flag_small_trees(lidR::remove_lasattribute(scene, "treeID")), "attribute 'treeID'")
  expect_error(flag_small_trees(lidR::remove_lasattribute(scene, "hag")), "attribute 'hag'")
  expect_error(flag_buffer(lidR::remove_lasattribute(scene, "treeID")), "attribute 'treeID'")
  expect_error(flag_buffer(scene, lidR::remove_lasattribute(seeds, "treeID")), "seeds must have an attribute 'treeID'")
  expect_error(flag_buffer(scene, seeds, 2), "'buffer' must be negative.")
  expect_error(flag_buffer(scene, seeds, "a"), "must be a number or an sf polygon")
  expect_error(flag_buffer(scene, seeds, multi), "simple polygon")
  expect_error(flag_buffer(scene, seeds, hole), "simple polygon")
  expect_error(flag_buffer(scene, seeds, other), "same CRS")
})

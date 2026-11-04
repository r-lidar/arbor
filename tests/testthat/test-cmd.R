test_that("check_args accepts valid command lines and rejects typos", {
  v <- arbor:::SEGMENT_VALUE_FLAGS
  b <- arbor:::SEGMENT_BOOL_FLAGS

  expect_true(arbor:::check_args(c("a.laz", "-buffer", "5", "-center", "-aoi", "x.shp"), v, b))
  expect_true(arbor:::check_args(c("a.laz", "-buffer", "-5"), v, b))
  expect_error(arbor:::check_args(c("a.laz", "-bound", "x.shp"), v, b), "Unknown option '-bound'")
  expect_error(arbor:::check_args(c("a.laz", "-buffer"), v, b), "requires a value")
  expect_error(arbor:::check_args(c("a.laz", "b.laz"), v, b), "Unexpected argument")
})

test_that("read_aoi accepts a simple polygon and rejects others", {
  skip_if_not_installed("sf")

  square <- function(x0) rbind(c(x0, 0), c(x0 + 1, 0), c(x0 + 1, 1), c(x0, 1), c(x0, 0))
  write_geom <- function(geom) {
    f <- tempfile(fileext = ".gpkg")
    suppressWarnings(sf::st_write(sf::st_sf(geometry = sf::st_sfc(geom, crs = 32631)), f, quiet = TRUE))
    f
  }

  f <- write_geom(sf::st_polygon(list(square(0))))
  expect_s3_class(arbor:::read_aoi(f), "sfc")

  f <- write_geom(sf::st_multipolygon(list(list(square(0)), list(square(5)))))
  expect_error(arbor:::read_aoi(f), "single polygon")

  hole <- sf::st_polygon(list(rbind(c(0, 0), c(4, 0), c(4, 4), c(0, 4), c(0, 0)),
                              rbind(c(1, 1), c(1, 2), c(2, 2), c(2, 1), c(1, 1))))
  expect_error(arbor:::read_aoi(write_geom(hole)), "single polygon")

  expect_error(arbor:::read_aoi("does_not_exist.shp"), "does not exist")
})

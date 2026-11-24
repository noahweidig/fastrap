test_that("input validation", {
  sq <- sf::st_sfc(sf::st_polygon(list(rbind(c(-100,40),c(-99.9,40),c(-99.9,40.1),c(-100,40)))), crs = 4326)
  expect_error(extract_rap(1985, sq), "1986")
  expect_error(extract_rap(2026, sq), "1986")
  expect_error(extract_rap(2020, sq, cover = "XYZ"), "cover")
  expect_error(extract_rap(2020, sf::st_sfc(sf::st_point(c(-100, 40)), crs = 4326)), "POLYGON")
})

test_that("tile keys cover bbox", {
  g <- sf::st_sfc(sf::st_buffer(sf::st_point(c(-100, 40)), 0.6), crs = 4326)
  bb <- fastrap:::feature_bboxes(g)
  expect_true(length(fastrap:::tiles_for(bb)) >= 4)
})

test_that("live extraction", {
  skip_on_cran()
  skip_if_offline("rangeland.ntsg.umt.edu")
  g <- sf::st_sfc(
    sf::st_buffer(sf::st_point(c(-107.5, 47.5)), 0.02),
    sf::st_buffer(sf::st_point(c(-107.4, 47.4)), 0.03),
    crs = 4326)
  cd <- tempfile()
  x <- extract_rap(2023, g, units = "ha", cache_dir = cd, quiet = TRUE)
  expect_equal(nrow(x), 12)
  expect_true(all(x$pct >= 0 & x$pct <= 100))
  w <- extract_rap(2023, g, cover = c("tre", "pfg"), wide = TRUE, cache_dir = cd, quiet = TRUE)
  expect_equal(names(w), c("id", "year", "pct_PFG", "area_PFG", "pct_TRE", "area_TRE"))
})

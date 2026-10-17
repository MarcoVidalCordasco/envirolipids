# Tests for npp_miami(). Expected values are computed directly from the
# published formulas, independently of the package code. The formulas give
# grams of dry matter, so they are multiplied by 0.45 to compare with the
# default output in grams of carbon.

test_that("npp_miami matches hand-computed values", {
  m <- npp_miami(mat = c(10, 25), map = c(500, 3000), id = c("a", "b"))
  expect_s3_class(m, "npp_climate")
  npp_t <- 3000 / (1 + exp(1.315 - 0.119 * c(10, 25)))
  npp_p <- 3000 * (1 - exp(-0.000664 * c(500, 3000)))
  expect_equal(m$npp, 0.45 * pmin(npp_t, npp_p))
  # 10 C / 500 mm: precipitation limits (npp_p ~ 848 vs npp_t ~ 1406).
  # 25 C / 3000 mm: temperature limits, narrowly (npp_t ~ 2521 vs npp_p ~ 2591).
  expect_equal(m$limiting, c("precipitation", "temperature"))
  expect_equal(m$id, c("a", "b"))
})

test_that("npp_miami can return dry matter instead of carbon", {
  mc <- npp_miami(mat = c(10, 25), map = c(500, 3000))
  md <- npp_miami(mat = c(10, 25), map = c(500, 3000), units = "dry_matter")
  expect_equal(mc$npp, 0.45 * md$npp)
  expect_equal(mc$npp_t, 0.45 * md$npp_t)
  expect_equal(mc$npp_p, 0.45 * md$npp_p)
  expect_equal(mc$limiting, md$limiting)
  expect_equal(attr(mc, "units"), "gC m-2 yr-1")
  expect_equal(attr(md, "units"), "g m-2 yr-1")
  expect_error(npp_miami(10, 500, units = "kg"), "should be one of")
})

test_that("npp_miami inputs are validated", {
  expect_error(npp_miami(1:3, 1:2), "same length")
  expect_error(npp_miami("10", 100), "numeric")
  expect_error(npp_miami(1:3, 1:3, id = 1:2), "`id`")
})

test_that("npp_miami results print with their units", {
  expect_output(print(npp_miami(10, 500)), "gC m-2 yr-1")
})

# Tests for the climate-based NPP models. Expected values are computed
# directly from the published formulas, independently of the package code.

test_that("npp_miami matches hand-computed values", {
  m <- npp_miami(mat = c(10, 25), map = c(500, 3000), id = c("a", "b"))
  expect_s3_class(m, "npp_climate")
  npp_t <- 3000 / (1 + exp(1.315 - 0.119 * c(10, 25)))
  npp_p <- 3000 * (1 - exp(-0.000664 * c(500, 3000)))
  expect_equal(m$npp, pmin(npp_t, npp_p))
  # 10 C / 500 mm: precipitation limits (npp_p ~ 848 vs npp_t ~ 1406).
  # 25 C / 3000 mm: temperature limits, narrowly (npp_t ~ 2521 vs npp_p ~ 2591).
  expect_equal(m$limiting, c("precipitation", "temperature"))
  expect_equal(m$id, c("a", "b"))
})

test_that("npp_thornthwaite matches hand-computed values", {
  # T = 10, P = 1000: L = 600, AET = 1000 / sqrt(0.9 + (1000 / 600)^2)
  # Here P > 0.316 L, so Turc's formula applies; AET is about 521 mm and NPP
  # about 1155 gC m-2 yr-1
  t <- npp_thornthwaite(mat = 10, map = 1000)
  aet <- 1000 / sqrt(0.9 + (1000 / 600)^2)
  expect_equal(t$aet, aet)
  expect_equal(t$npp, 3000 * (1 - exp(-0.0009695 * (aet - 20))))
})

test_that("npp_thornthwaite handles dry and cold climates", {
  # Sample 1: at 20 C, L = 1200, so P = 100 < 0.316 L and AET is capped at P.
  # Sample 2: at -30 C, L < 0 is floored at 0, so AET = 0 and NPP = 0.
  t <- npp_thornthwaite(mat = c(20, -30), map = c(100, 200))
  expect_equal(t$aet[1], 100) # P < 0.316 L, so AET = P
  expect_equal(t$npp[2], 0)
  expect_true(all(t$npp >= 0))
})

test_that("npp_thornthwaite gives 0, not NaN, with no rain in the cold", {
  # at -20 C, L is floored at 0; with P = 0 too, Turc's formula would be 0 / 0
  t <- npp_thornthwaite(mat = -20, map = 0)
  expect_equal(t$aet, 0)
  expect_equal(t$npp, 0)
})

test_that("climate inputs are validated", {
  expect_error(npp_miami(1:3, 1:2), "same length")
  expect_error(npp_thornthwaite("10", 100), "numeric")
  expect_error(npp_miami(1:3, 1:3, id = 1:2), "`id`")
})

test_that("climate autoplot returns a ggplot", {
  expect_s3_class(autoplot(npp_miami(1:5, 1:5 * 200)), "ggplot")
  expect_s3_class(autoplot(npp_thornthwaite(1:5, 1:5 * 200)), "ggplot")
})

# Tests for npp_thornthwaite(). Expected values are computed directly from the
# published formulas, independently of the package code. The formulas give
# grams of dry matter, so they are multiplied by 0.45 to compare with the
# default output in grams of carbon.

test_that("npp_thornthwaite matches hand-computed values", {
  # T = 10, P = 1000: L = 600, AET = 1000 / sqrt(0.9 + (1000 / 600)^2)
  # Here P > 0.316 L, so Turc's formula applies; AET is about 521 mm and NPP
  # about 1155 g m-2 yr-1 of dry matter
  t <- npp_thornthwaite(mat = 10, map = 1000)
  aet <- 1000 / sqrt(0.9 + (1000 / 600)^2)
  expect_equal(t$aet, aet)
  expect_equal(t$npp, 0.45 * 3000 * (1 - exp(-0.0009695 * (aet - 20))))
})

test_that("npp_thornthwaite can return dry matter instead of carbon", {
  tc <- npp_thornthwaite(mat = 10, map = 1000)
  td <- npp_thornthwaite(mat = 10, map = 1000, units = "dry_matter")
  expect_equal(tc$npp, 0.45 * td$npp)
  expect_equal(tc$aet, td$aet)
  expect_equal(attr(td, "units"), "g m-2 yr-1")
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

test_that("npp_thornthwaite inputs are validated", {
  expect_error(npp_thornthwaite(1:3, 1:2), "same length")
  expect_error(npp_thornthwaite("10", 100), "numeric")
  expect_error(npp_thornthwaite(1:3, 1:3, id = 1:2), "`id`")
})

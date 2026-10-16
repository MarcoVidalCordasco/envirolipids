# Tests for the input checks and transforms shared by all brGDGT functions.
# Internal objects such as BRGDGT_COLS are visible here because testthat runs
# the tests inside the package namespace.

test_that("fractions and percentages are both accepted", {
  # example_data is stored as fractions; scaling by 100 gives the same
  # samples as percentages, which must be converted back to identical fractions
  pct <- example_data
  pct[, BRGDGT_COLS] <- pct[, BRGDGT_COLS] * 100
  expect_equal(brgdgt_fractions(pct), brgdgt_fractions(example_data))
})

test_that("invalid brGDGT data is rejected", {
  # column 2 is one of the compounds, so dropping it must fail
  expect_error(brgdgt_fractions(example_data[, -2]), "missing brGDGT columns")
  # doubling one compound breaks the sum-to-1 check
  bad <- example_data
  bad$fIa <- bad$fIa * 2
  expect_error(brgdgt_fractions(bad), "must sum to 1")
  bad <- example_data
  bad$fIa[1] <- NA
  expect_error(brgdgt_fractions(bad), "missing values")
  expect_error(brgdgt_fractions(as.matrix(example_data)), "data frame")
})

# a defining property of the centred log-ratio transform
test_that("clr rows sum to zero", {
  clr <- brgdgt_clr(brgdgt_fractions(example_data))
  expect_equal(unname(rowSums(clr)), rep(0, nrow(example_data)))
})

test_that("methylation groups sum to one", {
  m <- brgdgt_methylation(example_data)
  expect_named(m, c("id", "tetra", "penta", "hexa"))
  expect_equal(m$id, example_data$ID)
  # drop the id column before summing the three groups
  expect_equal(rowSums(m[, -1]), rep(1, nrow(m)))
})

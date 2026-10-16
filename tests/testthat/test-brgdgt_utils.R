test_that("fractions and percentages are both accepted", {
  pct <- example_data
  pct[, BRGDGT_COLS] <- pct[, BRGDGT_COLS] * 100
  expect_equal(brgdgt_fractions(pct), brgdgt_fractions(example_data))
})

test_that("invalid brGDGT data is rejected", {
  expect_error(brgdgt_fractions(example_data[, -2]), "missing brGDGT columns")
  bad <- example_data
  bad$fIa <- bad$fIa * 2
  expect_error(brgdgt_fractions(bad), "must sum to 1")
  bad <- example_data
  bad$fIa[1] <- NA
  expect_error(brgdgt_fractions(bad), "missing values")
  expect_error(brgdgt_fractions(as.matrix(example_data)), "data frame")
})

test_that("clr rows sum to zero", {
  clr <- brgdgt_clr(brgdgt_fractions(example_data))
  expect_equal(unname(rowSums(clr)), rep(0, nrow(example_data)))
})

test_that("methylation groups sum to one", {
  m <- brgdgt_methylation(example_data)
  expect_named(m, c("id", "tetra", "penta", "hexa"))
  expect_equal(m$id, example_data$ID)
  expect_equal(rowSums(m[, -1]), rep(1, nrow(m)))
})

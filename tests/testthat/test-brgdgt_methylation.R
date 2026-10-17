# Tests for brgdgt_methylation()

test_that("methylation groups sum to one", {
  m <- brgdgt_methylation(example_data)
  expect_named(m, c("id", "tetra", "penta", "hexa"))
  expect_equal(m$id, example_data$ID)
  # drop the id column before summing the three groups
  expect_equal(rowSums(m[, -1]), rep(1, nrow(m)))
})

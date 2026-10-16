test_that("hanpp works with numeric vectors", {
  h <- hanpp(c(100, 200), c(300, 150))
  expect_s3_class(h, "hanpp")
  expect_equal(h$hanpp, c(200, -50))
  expect_equal(h$id, 1:2)
  expect_false("hanpp_lower" %in% names(h))
})

test_that("hanpp works with model objects and propagates intervals", {
  npp <- npp_brgdgt(example_data, "Lacustrine Sediment")
  climate <- npp_miami(mat = rep(12, nrow(example_data)),
                       map = rep(600, nrow(example_data)))
  h <- hanpp(npp, climate)
  expect_equal(h$id, example_data$ID)
  expect_equal(h$hanpp, climate$npp - npp$estimate)
  expect_equal(h$hanpp_lower, climate$npp - npp$upper)
  expect_equal(h$hanpp_upper, climate$npp - npp$lower)
  expect_equal(attr(h, "model"), "Miami")
  expect_s3_class(autoplot(h), "ggplot")
  expect_output(print(h), "Miami")
})

test_that("hanpp rejects mismatched or wrong inputs", {
  expect_error(hanpp(1:3, 1:2), "same number of samples")
  expect_error(hanpp("a", 1), "must be a numeric vector")
  mat <- temp_brgdgt(example_data, "Soil", model = load_npp_model())
  expect_error(hanpp(mat, seq_len(nrow(example_data))), "not NPP")
})

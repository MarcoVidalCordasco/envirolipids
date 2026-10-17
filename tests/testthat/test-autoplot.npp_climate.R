# Tests for autoplot() on npp_miami() and npp_thornthwaite() results

test_that("climate autoplot returns a ggplot", {
  expect_s3_class(autoplot(npp_miami(1:5, 1:5 * 200)), "ggplot")
  expect_s3_class(autoplot(npp_thornthwaite(1:5, 1:5 * 200)), "ggplot")
})

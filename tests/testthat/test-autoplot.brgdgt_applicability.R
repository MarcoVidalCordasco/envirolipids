# Tests for autoplot() on npp_brgdgt_applicable() results

check <- npp_brgdgt_applicable(example_data, n_neighbors = 3)

test_that("all plot types return a ggplot", {
  for (type in c("mahalanobis", "pca", "ternary")) {
    expect_s3_class(autoplot(check, type = type), "ggplot")
  }
  # the Mahalanobis plot can also be drawn against a numeric x, such as age
  p <- autoplot(check, type = "mahalanobis", x = example_data$Age)
  expect_true(is.numeric(p$data$x))
  # the map needs the suggested maps package
  skip_if_not_installed("maps")
  expect_s3_class(autoplot(check, type = "map"), "ggplot")
})

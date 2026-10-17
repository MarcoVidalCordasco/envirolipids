# Tests for autoplot() on npp_brgdgt() and temp_brgdgt() predictions

test_that("prediction autoplot returns a ggplot", {
  npp <- npp_brgdgt(example_data, "Soil")
  expect_s3_class(autoplot(npp), "ggplot")
  expect_s3_class(autoplot(npp, x = example_data$Age), "ggplot")
  # x must align with the samples
  expect_error(autoplot(npp, x = 1:3), "one value per sample")
})

test_that("autoplot draws samples by default and a series along numeric x", {
  npp <- npp_brgdgt(modern_example_data, modern_example_data$Sampletype)
  # by sample: ids on a discrete axis, error bars, no connecting line
  by_sample <- autoplot(npp)
  geoms <- vapply(by_sample$layers, function(l) class(l$geom)[1], "")
  expect_true(is.factor(by_sample$data$x))
  expect_equal(levels(by_sample$data$x), modern_example_data$ID)
  expect_true("GeomErrorbar" %in% geoms)
  expect_false("GeomLine" %in% geoms)
  # along a numeric x: a line with a ribbon
  series <- autoplot(npp, x = modern_example_data$Latitude)
  geoms <- vapply(series$layers, function(l) class(l$geom)[1], "")
  expect_true(all(c("GeomRibbon", "GeomLine") %in% geoms))
})

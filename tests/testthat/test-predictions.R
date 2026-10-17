# Tests for npp_brgdgt() and temp_brgdgt(). The regression check against the
# pre-refactor code (identical estimates and intervals) was done once during
# the refactor; these tests cover behaviour that must keep holding.

test_that("npp_brgdgt returns a brgdgt_prediction", {
  npp <- npp_brgdgt(example_data, "Lacustrine Sediment")
  expect_s3_class(npp, "brgdgt_prediction")
  expect_named(npp, c("id", "sample_type", "estimate", "lower", "upper"))
  expect_equal(nrow(npp), nrow(example_data))
  expect_equal(npp$id, example_data$ID)
  # The forest mean need not lie between the 2.5% and 97.5% tree quantiles in
  # general, but for a well-behaved forest it should
  expect_true(all(npp$lower <= npp$estimate & npp$estimate <= npp$upper))
  expect_equal(attr(npp, "response"), "NPP")
})

# The interval only adds columns; it must not change the point estimate
test_that("npp_brgdgt without ci gives the same estimates", {
  with_ci <- npp_brgdgt(example_data, "Soil")
  without <- npp_brgdgt(example_data, "Soil", ci = FALSE)
  expect_named(without, c("id", "sample_type", "estimate"))
  expect_equal(without$estimate, with_ci$estimate)
})

test_that("brGDGT columns are matched by name, not position", {
  # reversing the column order would break any code selecting by position
  shuffled <- example_data[, rev(names(example_data))]
  expect_equal(npp_brgdgt(shuffled, "Peat")$estimate,
               npp_brgdgt(example_data, "Peat")$estimate)
})

test_that("sample_type is validated and can vary by sample", {
  expect_error(npp_brgdgt(example_data, "Ocean"), "Unknown `sample_type`")
  expect_error(npp_brgdgt(example_data, c("Soil", "Peat")), "length 1")
  # Alternate sample types: each sample must get the same prediction as if
  # its type had been given for the whole dataset
  types <- rep(c("Soil", "Peat"), length.out = nrow(example_data))
  mixed <- npp_brgdgt(example_data, types)
  expect_equal(mixed$estimate[types == "Peat"],
               npp_brgdgt(example_data, "Peat")$estimate[types == "Peat"])
})

# These types are pooled into one class before reaching the model
test_that("aquatic SPM types share a prediction", {
  expect_equal(npp_brgdgt(example_data, "Lacustrine SPM")$estimate,
               npp_brgdgt(example_data, "Riverine Sediment and SPM")$estimate)
})

test_that("ids default to row numbers", {
  no_id <- example_data[, names(example_data) != "ID"]
  expect_equal(npp_brgdgt(no_id, "Soil")$id, seq_len(nrow(no_id)))
})

test_that("temp_brgdgt errors clearly without a model", {
  # once the temperature model ships, this test no longer applies
  skip_if(nzchar(system.file("extdata", "model_brGDGT-MAT.rds",
                             package = "envirolipids")))
  expect_error(temp_brgdgt(example_data, "Soil"),
               "temperature model is not yet included")
})

test_that("temp_brgdgt works with a supplied model", {
  # No temperature model exists yet, so borrow the NPP model to check the
  # plumbing: the values are meaningless, but the labelling must be right
  mat <- temp_brgdgt(example_data, "Soil", model = load_npp_model())
  expect_s3_class(mat, "brgdgt_prediction")
  expect_equal(attr(mat, "response"), "MAT")
})

test_that("prediction autoplot returns a ggplot", {
  npp <- npp_brgdgt(example_data, "Soil")
  expect_s3_class(autoplot(npp), "ggplot")
  expect_s3_class(autoplot(npp, x = example_data$Age), "ggplot")
  # x must align with the samples
  expect_error(autoplot(npp, x = 1:3), "one value per sample")
  expect_output(print(npp), "brGDGT-based NPP predictions")
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

test_that("npp_brgdgt returns a brgdgt_prediction", {
  npp <- npp_brgdgt(example_data, "Lacustrine Sediment")
  expect_s3_class(npp, "brgdgt_prediction")
  expect_named(npp, c("id", "sample_type", "estimate", "lower", "upper"))
  expect_equal(nrow(npp), nrow(example_data))
  expect_equal(npp$id, example_data$ID)
  expect_true(all(npp$lower <= npp$estimate & npp$estimate <= npp$upper))
  expect_equal(attr(npp, "response"), "NPP")
})

test_that("npp_brgdgt without ci gives the same estimates", {
  with_ci <- npp_brgdgt(example_data, "Soil")
  without <- npp_brgdgt(example_data, "Soil", ci = FALSE)
  expect_named(without, c("id", "sample_type", "estimate"))
  expect_equal(without$estimate, with_ci$estimate)
})

test_that("brGDGT columns are matched by name, not position", {
  shuffled <- example_data[, rev(names(example_data))]
  expect_equal(npp_brgdgt(shuffled, "Peat")$estimate,
               npp_brgdgt(example_data, "Peat")$estimate)
})

test_that("sample_type is validated and can vary by sample", {
  expect_error(npp_brgdgt(example_data, "Ocean"), "Unknown `sample_type`")
  expect_error(npp_brgdgt(example_data, c("Soil", "Peat")), "length 1")
  types <- rep(c("Soil", "Peat"), length.out = nrow(example_data))
  mixed <- npp_brgdgt(example_data, types)
  expect_equal(mixed$estimate[types == "Peat"],
               npp_brgdgt(example_data, "Peat")$estimate[types == "Peat"])
})

test_that("aquatic SPM types share a prediction", {
  expect_equal(npp_brgdgt(example_data, "Lacustrine SPM")$estimate,
               npp_brgdgt(example_data, "Riverine Sediment and SPM")$estimate)
})

test_that("ids default to row numbers", {
  no_id <- example_data[, names(example_data) != "ID"]
  expect_equal(npp_brgdgt(no_id, "Soil")$id, seq_len(nrow(no_id)))
})

test_that("temp_brgdgt errors clearly without a model", {
  skip_if(nzchar(system.file("extdata", "model_brGDGT-MAT.rds",
                             package = "envirolipids")))
  expect_error(temp_brgdgt(example_data, "Soil"),
               "temperature model is not yet included")
})

test_that("temp_brgdgt works with a supplied model", {
  mat <- temp_brgdgt(example_data, "Soil", model = load_npp_model())
  expect_s3_class(mat, "brgdgt_prediction")
  expect_equal(attr(mat, "response"), "MAT")
})

test_that("prediction autoplot returns a ggplot", {
  npp <- npp_brgdgt(example_data, "Soil")
  expect_s3_class(autoplot(npp), "ggplot")
  expect_s3_class(autoplot(npp, x = example_data$Age), "ggplot")
  expect_error(autoplot(npp, x = 1:3), "one value per sample")
  expect_output(print(npp), "brGDGT-based NPP predictions")
})

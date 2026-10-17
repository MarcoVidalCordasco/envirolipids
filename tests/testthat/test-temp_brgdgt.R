# Tests for temp_brgdgt()

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

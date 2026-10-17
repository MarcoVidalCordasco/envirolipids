# Tests for npp_brgdgt_applicable() and its print/autoplot methods

# Run the (relatively slow) check once and share it between the tests below
check <- npp_brgdgt_applicable(example_data, n_neighbors = 3)

test_that("npp_brgdgt_applicable returns per-sample diagnostics", {
  expect_s3_class(check, "brgdgt_applicability")
  s <- as.data.frame(check)
  expect_named(s, c("id", "mahalanobis", "p_value", "within_domain",
                    "nn_distance"))
  expect_equal(s$id, example_data$ID)
  # the domain flag must agree with the reported cut-off
  expect_equal(s$within_domain, s$mahalanobis <= check$threshold_value)
  expect_true(all(s$p_value >= 0 & s$p_value <= 1))
  expect_equal(check$method, "pca")
})

test_that("the cut-off keeps the requested share of training samples", {
  # Training samples are, by definition, in the domain: with the empirical
  # cut-off, the proportion of them inside it matches `threshold`
  reference <- load_reference_data()
  ref_data <- reference$fractions
  ref_data$ID <- seq_len(nrow(ref_data))
  ref_check <- npp_brgdgt_applicable(ref_data, threshold = 0.9)
  expect_equal(mean(ref_check$samples$within_domain), 0.9, tolerance = 0.01)
  # a sample on the cut-off has 10% of training samples at least as far out
  expect_equal(mean(ref_check$samples$p_value <= 0.1), 0.1, tolerance = 0.01)
})

test_that("neighbours are sorted and consistent with nn_distance", {
  # 3 neighbours per sample, the first of which is the nearest overall
  nb <- check$neighbours
  expect_equal(nrow(nb), 3 * nrow(example_data))
  expect_equal(nb$distance[nb$rank == 1], check$samples$nn_distance)
  # neighbours of a sample are listed nearest first
  expect_true(all(diff(nb$distance[nb$id == 1]) >= 0))
})

test_that("regularized method works", {
  reg <- npp_brgdgt_applicable(example_data, use_pca = FALSE)
  expect_equal(reg$method, "regularized")
  # without PCA, the degrees of freedom are all 15 CLR components
  expect_equal(reg$n_pcs, length(BRGDGT_COLS))
})

test_that("results for a sample do not depend on the rest of the batch", {
  # The PCA is fitted on the training data only, so checking one sample on its
  # own must give exactly what it got as part of the full dataset. Also
  # exercises the single-row edge case.
  one <- npp_brgdgt_applicable(example_data[1, ])
  expect_equal(nrow(one$samples), 1)
  expect_equal(one$samples, check$samples[1, ], ignore_attr = TRUE)
  expect_equal(one$n_pcs, check$n_pcs)
})

test_that("print and autoplot work", {
  expect_output(print(check), "Samples within domain: ")
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

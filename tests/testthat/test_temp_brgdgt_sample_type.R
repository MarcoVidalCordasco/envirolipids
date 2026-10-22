# Tests for published brGDGT temperature calibrations.
# Positive abundances ensure every test index is defined.

temp_test_data <- function() {
  x <- c(fIa = 0.30, fIb = 0.06, fIc = 0.02,
         fIIa = 0.20, fIIa_ = 0.08,
         fIIb = 0.05, fIIb_ = 0.03,
         fIIc = 0.02, fIIc_ = 0.01,
         fIIIa = 0.11, fIIIa_ = 0.04,
         fIIIb = 0.03, fIIIb_ = 0.02,
         fIIIc = 0.02, fIIIc_ = 0.01)
  x <- x / sum(x)
  data.frame(ID = "test1", as.list(x), check.names = FALSE)
}

test_that("all advertised temperature methods execute", {
  x <- temp_test_data()
  catalog <- temp_brgdgt_methods()
  expect_equal(nrow(catalog), 18L)
  expect_equal(length(unique(catalog$method)), 18L)
  for (m in catalog$method) {
    sample_type <- catalog$sample_group[catalog$method == m]
    y <- suppressWarnings(temp_brgdgt(x, sample_type = sample_type,
                                      method = m))
    expect_s3_class(y, "brgdgt_prediction")
    expect_equal(nrow(y), 1L)
    expect_equal(y$method, m)
    expect_true(y$target %in% c("MAT", "MAF", "GDD0"))
    expect_identical(attr(y, "response"), y$target[1L])
  }
})

test_that("Naafs, De Jonge and Russell equations match their coefficients", {
  x <- temp_test_data()
  f <- brgdgt_fractions(x)
  mbt <- (f$fIa + f$fIb + f$fIc) /
    (f$fIa + f$fIb + f$fIc + f$fIIa + f$fIIb + f$fIIc + f$fIIIa)
  get <- function(method) {
    catalogue <- temp_brgdgt_methods()
    sample_type <- catalogue$sample_group[catalogue$method == method]
    temp_brgdgt(x, sample_type = sample_type, method = method)$estimate
  }
  
  expect_equal(get("dejonge2014_mbt"), -8.57 + 31.45 * mbt)
  expect_equal(get("dejonge2014_mlr"),
               7.17 + 17.1 * f$fIa + 25.9 * f$fIb +
                 34.4 * f$fIc - 28.6 * f$fIIa)
  expect_equal(get("naafs2017_global"), 40.01 * mbt - 15.25)
  expect_equal(get("naafs2017_gdd"), 14344.3 * mbt - 4997.5)
  expect_equal(get("naafs2017_mlr"),
               19.8 * f$fIa + 31.1 * f$fIb - 23.4 * f$fIIa + 4.32)
  expect_equal(get("naafs2017_gdd_mlr"),
               6152.9 * f$fIa + 8272.1 * f$fIb -
                 8015.8 * f$fIIa + 2509.4)
  expect_equal(get("russell2018_mbt"), -1.21 + 32.42 * mbt)
  expect_equal(get("raberg2021_mbt_maf"), -0.5 + 30.4 * mbt)
  expect_equal(get("zhu2025_1"), 28.830 * mbt - 0.268)
  expect_equal(get("zhu2025_2"),
               19.901 * f$fIa + 43.622 * f$fIb -
                 16.041 * f$fIIa + 9.706 * f$fIIa_ -
                 7.293 * f$fIIIa - 5.578 * f$fIIIa_ + 8.233)
  expect_equal(get("zhu2025_6"), 23.966 * f$fIa + 70.362 * f$fIb + 2.166)
})

test_that("original Peterse 2012 equation uses MBT prime and log10 CBT", {
  x <- temp_test_data()
  f <- brgdgt_fractions(x)
  mbtp <- (f$fIa + f$fIb + f$fIc) /
    (f$fIa + f$fIb + f$fIc + f$fIIa + f$fIIa_ +
       f$fIIb + f$fIIb_ + f$fIIc + f$fIIc_ + f$fIIIa + f$fIIIa_)
  cbt <- -log10((f$fIb + f$fIIb + f$fIIb_) /
                  (f$fIa + f$fIIa + f$fIIa_))
  expect_equal(temp_brgdgt(x, sample_type = "soils",
                           method = "peterse2012")$estimate,
               0.81 - 5.67 * cbt + 31.0 * mbtp)
})

test_that("Naafs restricted calibration refuses high IR6me", {
  x <- temp_test_data()
  # Force strong six-methyl dominance with an otherwise valid composition.
  cols6 <- c("fIIa_", "fIIb_", "fIIc_", "fIIIa_", "fIIIb_", "fIIIc_")
  x[cols6] <- x[cols6] * 15
  cols <- BRGDGT_COLS
  x[cols] <- x[cols] / sum(unlist(x[cols]))
  res <- suppressWarnings(temp_brgdgt(x, sample_type = "mineral_soils",
                                      method = "naafs2017_restricted"))
  expect_true(res$IR6me >= 0.5)
  expect_false(res$applicable)
  expect_true(is.na(res$estimate))
})

test_that("MAF is not mislabeled MAT or degree days", {
  x <- temp_test_data()
  maf <- temp_brgdgt(x, sample_type = "lake_sediments",
                     method = "raberg2021_eq10")
  expect_equal(maf$target, "MAF")
  expect_equal(maf$calibrated_on, "Lake sediment")
  gdd <- temp_brgdgt(x, sample_type = "mineral_soils",
                     method = "naafs2017_gdd_mlr_restricted")
  expect_equal(gdd$target, "GDD0")
  expect_equal(gdd$units, "degC-days/year")
})

test_that("sample type filters the available methods", {
  expect_equal(nrow(temp_brgdgt_methods("soils")), 3L)
  expect_equal(nrow(temp_brgdgt_methods("mineral_soils")), 11L)
  expect_equal(nrow(temp_brgdgt_methods("lake_sediments")), 7L)
  expect_equal(nrow(temp_brgdgt_methods()), 18L)
  
  expect_true("naafs2017_global" %in%
                temp_brgdgt_methods("mineral_soils")$method)
  expect_false("naafs2017_global" %in%
                 temp_brgdgt_methods("soils")$method)
  expect_true("dejonge2014_mbt" %in%
                temp_brgdgt_methods("mineral_soils")$method)
  expect_false("zhu2025_1" %in%
                 temp_brgdgt_methods("mineral_soils")$method)
})

test_that("methods incompatible with the selected material fail", {
  x <- temp_test_data()
  expect_error(temp_brgdgt(x, sample_type = "soils",
                           method = "zhu2025_1"),
               "not available")
  expect_error(temp_brgdgt(x, sample_type = "lake_sediments",
                           method = "naafs2017_global"),
               "not available")
  expect_error(temp_brgdgt(x, sample_type = "soils",
                           method = "naafs2017_global"),
               "not available")
  expect_error(temp_brgdgt(x, sample_type = "mineral_soils"),
               "Please choose")
  expect_error(temp_brgdgt_methods("peat"), "Unknown")
})

test_that("friendly material labels and prediction metadata work", {
  x <- temp_test_data()
  out <- temp_brgdgt(x, sample_type = "mineral soils",
                     method = "peterse2012")
  expect_equal(out$sample_type, "mineral_soils")
  expect_true(all(out$environment_match))
  expect_equal(out$calibrated_on, "Soils")
  lake <- temp_brgdgt(x, sample_type = "lake sediment",
                      method = "zhu2025_1")
  expect_equal(lake$target, "MAF")
  expect_equal(lake$sample_type, "lake_sediments")
})

test_that("legacy random forest method accepts a supplied fitted model", {
  x <- temp_test_data()
  fit <- load_npp_model()  # exercise plumbing; not a real MAT model
  out <- temp_brgdgt(x, sample_type = "soils",
                     model = fit, method = "random_forest", ci = FALSE)
  expect_s3_class(out, "brgdgt_prediction")
  expect_equal(attr(out, "response"), "MAT")
})

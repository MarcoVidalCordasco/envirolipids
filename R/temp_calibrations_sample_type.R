# Published brGDGT temperature calibrations.
# Targets: MAT (mean annual AIR temperature, degC), MAF (mean AIR
# temperature of months above 0 degC, degC), GDD0 (degree-days/year).
# Reference material means what was sampled to FIT the calibration.

# Safe, row-wise ratio: undefined denominators become NA, not Inf or NaN.
.temp_ratio <- function(x, denominator) {
  out <- x / denominator
  out[!is.finite(out)] <- NA_real_
  out
}

# Named 15-compound fractions are passed as f, with rows normalised to 1.
.temp_indices <- function(f) {
  v <- function(cols) rowSums(f[, cols, drop = FALSE])
  tetra <- v(c("fIa", "fIb", "fIc"))
  five_penta_hexa <- v(c("fIIa", "fIIb", "fIIc",
                         "fIIIa", "fIIIb", "fIIIc"))
  six_penta_hexa <- v(c("fIIa_", "fIIb_", "fIIc_",
                        "fIIIa_", "fIIIb_", "fIIIc_"))
  mbt5me <- .temp_ratio(tetra,
                        v(c("fIa", "fIb", "fIc",
                            "fIIa", "fIIb", "fIIc", "fIIIa")))
  # Original MBT' (Peterse et al., 2012): includes 6-methyl isomers,
  # but excludes IIIb/c, including their 6-methyl isomers.
  mbtprime <- .temp_ratio(tetra,
                          v(c("fIa", "fIb", "fIc", "fIIa", "fIIa_",
                              "fIIb", "fIIb_", "fIIc", "fIIc_",
                              "fIIIa", "fIIIa_")))
  # Original CBT is a BASE-10 log (not natural log).
  cbt_fraction <- .temp_ratio(f$fIb + f$fIIb + f$fIIb_,
                              f$fIa + f$fIIa + f$fIIa_)
  cbt <- -log10(cbt_fraction)
  cbt[!is.finite(cbt)] <- NA_real_
  ir6me <- .temp_ratio(six_penta_hexa,
                       five_penta_hexa + six_penta_hexa)
  list(mbt5me = mbt5me, mbtprime = mbtprime,
       cbt = cbt, ir6me = ir6me)
}

# Fractions within the structural methylation groups of Raberg et al.
# (2021), Table A1. They differ from conventional 15-compound fractions.
.temp_raberg_meth <- function(f) {
  a <- f$fIa + f$fIIa + f$fIIIa + f$fIIa_ + f$fIIIa_
  b <- f$fIb + f$fIIb + f$fIIIb + f$fIIb_ + f$fIIIb_
  c <- f$fIc + f$fIIc + f$fIIIc + f$fIIc_ + f$fIIIc_
  list(ib = .temp_ratio(f$fIb, b),
       iia = .temp_ratio(f$fIIa, a),
       iib = .temp_ratio(f$fIIb, b),
       iic = .temp_ratio(f$fIIc, c),
       iiia = .temp_ratio(f$fIIIa, a),
       iiib = .temp_ratio(f$fIIIb, b))
}

# Each row is an IMPLEMENTED method only. RMSEs use the original paper's
# definition; do not interpret as per-sample 95% confidence intervals.
.temp_catalogue <- function() {
  data.frame(
    method = c(
      "peterse2012", "dejonge2014_mbt", "dejonge2014_mlr",
      "naafs2017_global", "naafs2017_gdd", "naafs2017_mlr",
      "naafs2017_gdd_mlr", "naafs2017_restricted",
      "naafs2017_gdd_restricted", "naafs2017_mlr_restricted",
      "naafs2017_gdd_mlr_restricted",
      "russell2018_mbt", "raberg2021_mbt_maf",
      "raberg2021_eq10", "raberg2021_eq11",
      "zhu2025_1", "zhu2025_2", "zhu2025_6"),
    target = c(rep("MAT", 4), "GDD0", "MAT", "GDD0", "MAT", "GDD0",
               "MAT", "GDD0", "MAT", rep("MAF", 6)),
    units = c(rep("degC", 4), "degC-days/year", "degC", "degC-days/year",
              "degC", "degC-days/year", "degC", "degC-days/year",
              rep("degC", 7)),
    # Published calibration material (precise category reported to users).
    calibrated_on = c(rep("Soils", 3), rep("Mineral soils", 8),
                      rep("Lake sediment", 4),
                      "Freshwater lake sediment", "Freshwater lake sediment",
                      "Freshwater and non-freshwater lake sediment"),
    # Used only to decide which methods are offered for a sample type.
    # Mineral-soil users may also use the earlier general soil equations.
    sample_group = c(rep("soils", 3), rep("mineral_soils", 8),
                     rep("lake_sediments", 7)),
    publication = c("Peterse et al. (2012)",
                    rep("De Jonge et al. (2014)", 2),
                    rep("Naafs et al. (2017)", 8),
                    "Russell et al. (2018)",
                    rep("Raberg et al. (2021)", 3),
                    rep("Zhu et al. (2025)", 3)),
    equation = c("Eq. 3 in Naafs (2017)", "Eq. 5 in Naafs (2017)",
                 "Eq. 6 in Naafs (2017)",
                 "Eq. 9", "Eq. 10", "Eq. 11", "Eq. 12",
                 "Eq. 14", "Eq. 15", "Eq. 16", "Eq. 17",
                 "MBT'5me", "MBT'5me", "Eq. 10", "Eq. 11",
                 "Calibration 1", "Calibration 2", "Calibration 6"),
    rmse = c(5.0, 4.9, 4.7, 5.3, 1779, 4.7, 1319,
             4.1, 1326, 3.8, 1079, 2.14, NA_real_, 2.14, 1.97,
             2.65, 2.33, NA_real_),
    rmse_context = c(rep("calibration", 11), "calibration", NA_character_,
                     "calibration", "calibration", "test set", "test set",
                     "2.87 freshwater; 3.62 non-freshwater"),
    requires_ir6me_below_0.5 = c(rep(FALSE, 7), TRUE, TRUE, TRUE, TRUE,
                                 rep(FALSE, 7)),
    notes = c(
      "Original 2012 soil MBT'/CBT formulation, via Naafs (2017)",
      "Original study, via Naafs (2017)",
      "Original study, via Naafs (2017)",
      rep("Soil calibration; GDD0 is not a mean temperature", 8),
      "Lacustrine MAT from East African lake sediments (n=65)",
      "Index-only MAF relationship in Raberg discussion; not their recommended Eq. 10",
      "Raberg methylation-group fractional-abundance quadratic regression",
      "Raberg full 15-compound fractional-abundance quadratic regression",
      "Freshwater lake sediment; MAF, not MAT",
      "Freshwater lake sediment; MAF, not MAT",
      "Fitted on mixed lake dataset; performance differs by salinity"),
    stringsAsFactors = FALSE
  )
}

# Equation dispatcher: internal implementation is separate from UI/metadata.
.predict_temp_equation <- function(f, method, idx) {
  a <- f$fIa; b <- f$fIb; c <- f$fIc
  iia <- f$fIIa; iiia <- f$fIIIa
  iia6 <- f$fIIa_; iiia6 <- f$fIIIa_
  m <- idx$mbt5me
  if (method == "peterse2012") {
    return(0.81 - 5.67 * idx$cbt + 31.0 * idx$mbtprime)
  }
  if (method == "dejonge2014_mbt") return(-8.57 + 31.45 * m)
  if (method == "dejonge2014_mlr") {
    return(7.17 + 17.1 * a + 25.9 * b + 34.4 * c - 28.6 * iia)
  }
  if (method == "naafs2017_global") return(40.01 * m - 15.25)
  if (method == "naafs2017_gdd") return(14344.3 * m - 4997.5)
  if (method == "naafs2017_mlr") {
    return(19.8 * a + 31.1 * b - 23.4 * iia + 4.32)
  }
  if (method == "naafs2017_gdd_mlr") {
    return(6152.9 * a + 8272.1 * b - 8015.8 * iia + 2509.4)
  }
  if (method == "naafs2017_restricted") return(39.09 * m - 14.50)
  if (method == "naafs2017_gdd_restricted") {
    return(13498.8 * m - 4444.5)
  }
  if (method == "naafs2017_mlr_restricted") {
    return(14.7 * a - 31.7 * iia + 10.0)
  }
  if (method == "naafs2017_gdd_mlr_restricted") {
    return(4881 * a - 10112 * iia + 3942)
  }
  if (method == "russell2018_mbt") return(-1.21 + 32.42 * m)
  if (method == "raberg2021_mbt_maf") return(-0.5 + 30.4 * m)
  if (method == "raberg2021_eq10") {
    r <- .temp_raberg_meth(f)
    return(92.9 + 63.84 * r$ib^2 - 130.51 * r$ib -
             28.77 * r$iia^2 - 72.28 * r$iib^2 - 5.88 * r$iic^2 +
             20.89 * r$iiia^2 - 40.54 * r$iiia - 80.47 * r$iiib)
  }
  if (method == "raberg2021_eq11") {
    return(-8.06 + 37.52 * a - 266.83 * b^2 + 133.42 * b +
             100.85 * iia6^2 + 58.15 * iiia6^2 + 12.79 * iiia)
  }
  if (method == "zhu2025_1") return(28.830 * m - 0.268)
  if (method == "zhu2025_2") {
    return(19.901 * a + 43.622 * b - 16.041 * iia +
             9.706 * iia6 - 7.293 * iiia - 5.578 * iiia6 + 8.233)
  }
  if (method == "zhu2025_6") return(23.966 * a + 70.362 * b + 2.166)
  stop("Unimplemented temperature calibration: ", method, call. = FALSE)
}

# Apply and return one temperature metric per sample.
predict_temp_calibration <- function(fractions, method) {
  meta <- .temp_catalogue()
  match_row <- match(method, meta$method)
  if (is.na(match_row)) stop("Unknown temperature method: ", method, call. = FALSE)
  ref <- meta[match_row, , drop = FALSE]
  # The input helper allows a 1% rounding tolerance, so correct that here
  # before calculating individual compound FAs in regression equations.
  f <- as.data.frame(sweep(as.matrix(fractions), 1L,
                           rowSums(fractions), "/"))
  idx <- .temp_indices(f)
  estimate <- .predict_temp_equation(f, method, idx)
  applicable <- if (ref$requires_ir6me_below_0.5) {
    !is.na(idx$ir6me) & idx$ir6me < 0.5
  } else rep(TRUE, nrow(f))
  estimate[!applicable] <- NA_real_
  data.frame(estimate = estimate, target = ref$target,
             units = ref$units, calibrated_on = ref$calibrated_on,
             calibration_group = ref$sample_group,
             method = method, publication = ref$publication,
             equation = ref$equation, MBT5me = idx$mbt5me,
             IR6me = idx$ir6me, applicable = applicable,
             rmse = ref$rmse, stringsAsFactors = FALSE)
}

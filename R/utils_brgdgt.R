# Internal helpers for reading and transforming brGDGT data.

# Extract the brGDGT columns as fractional abundances, accepting either
# fractions (rows sum to 1) or percentages (rows sum to 100)
brgdgt_fractions <- function(data) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame", call. = FALSE)
  }
  # all 15 compounds are needed, as the models use every one of them
  missing_cols <- setdiff(BRGDGT_COLS, names(data))
  if (length(missing_cols) > 0) {
    stop("`data` is missing brGDGT columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  # as.data.frame() drops tibble behaviour, so that `[` and rowSums() below
  # behave the same whatever kind of data frame the user passed in
  x <- as.data.frame(data)[, BRGDGT_COLS]
  if (!all(vapply(x, is.numeric, logical(1)))) {
    stop("brGDGT columns must be numeric", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("brGDGT columns must not contain missing values", call. = FALSE)
  }
  if (any(x < 0)) {
    stop("brGDGT abundances must not be negative", call. = FALSE)
  }
  # Work out the units from the row totals. The tolerances (+/- 1% of the
  # total) allow for rounding in published tables; the training data
  # themselves sum to between 99.8 and 100.2. We do not rescale rows to sum
  # exactly to 1, so that new data are treated the same way as the training
  # data were.
  totals <- rowSums(x)
  if (all(abs(totals - 100) <= 1)) {
    x <- x / 100
  } else if (!all(abs(totals - 1) <= 0.01)) {
    stop("brGDGT abundances must sum to 1 (fractions) or 100 (percentages) ",
         "for every sample", call. = FALSE)
  }
  x
}

# Centred log-ratio transform of fractional abundances, with a small offset
# for zeros (as used to train the models)
brgdgt_clr <- function(fractions) {
  # clr(x) = log(x) - mean(log(x)), computed row by row. Zeros have no
  # logarithm, so a tiny offset (1e-6, i.e. 0.0001%) is added to every
  # compound first. The offset must stay the same as the one used to build the
  # training data (data-raw/make_reference_data.R), and it is only consistent
  # because inputs are always fractions at this point.
  log_comp <- log(fractions + 1e-6)
  # rowMeans() returns one value per sample; subtracting it from the matrix
  # recycles it down the columns, i.e. it centres each row
  as.data.frame(log_comp - rowMeans(log_comp))
}

# Sample identifiers: the user's id column if present, else row numbers
get_ids <- function(data, id_col) {
  if (!is.null(id_col) && id_col %in% names(data)) {
    data[[id_col]]
  } else {
    seq_len(nrow(data))
  }
}

# Internal workhorse for brgdgt_methylation(), also used on the training data
# (which has no id column) by npp_brgdgt_applicable()
methylation_groups <- function(fractions) {
  # one column per group, holding the summed abundance of its compounds
  groups <- vapply(METHYLATION_GROUPS,
                   function(cols) rowSums(fractions[, cols, drop = FALSE]),
                   numeric(nrow(fractions)))
  # vapply() simplifies a single sample to a plain vector, so rebuild the
  # matrix explicitly to always get one row per sample
  groups <- matrix(groups, nrow = nrow(fractions),
                   dimnames = list(NULL, names(METHYLATION_GROUPS)))
  # rescale so the three groups sum exactly to 1, as needed for a ternary plot
  as.data.frame(groups / rowSums(groups))
}



# MBT'5me index (Naafs et al., 2017, Eq. 4) This is used to test whether some specific predictive
# formulas of MAT can be used:

brgdgt_mbt5me <- function(fractions) {
  
  numerator <- rowSums(
    fractions[, c("fIa", "fIb", "fIc"), drop = FALSE]
  )
  
  denominator <- rowSums(
    fractions[, c(
      "fIa", "fIb", "fIc",
      "fIIa", "fIIb", "fIIc", "fIIIa"
    ), drop = FALSE]
  )
  
  result <- numerator / denominator
  result[denominator == 0] <- NA_real_
  
  result
}

# IR6me: relative abundance of 6-methyl isomers
# among penta- and hexamethylated brGDGTs
brgdgt_ir6me <- function(fractions) {
  
  six_methyl <- rowSums(
    fractions[, c(
      "fIIa_", "fIIb_", "fIIc_",
      "fIIIa_", "fIIIb_", "fIIIc_"
    ), drop = FALSE]
  )
  
  five_methyl <- rowSums(
    fractions[, c(
      "fIIa", "fIIb", "fIIc",
      "fIIIa", "fIIIb", "fIIIc"
    ), drop = FALSE]
  )
  
  denominator <- five_methyl + six_methyl
  
  result <- six_methyl / denominator
  result[denominator == 0] <- NA_real_
  
  result
}

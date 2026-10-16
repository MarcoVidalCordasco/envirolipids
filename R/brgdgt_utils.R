#' @importFrom ggplot2 autoplot .data
#' @importFrom randomForest randomForest
NULL

# Re-export autoplot() so that users can call autoplot() on our objects
# without having to attach ggplot2 themselves
#' @export
ggplot2::autoplot

# The 15 brGDGT compounds used by the models. All functions select the
# compounds from user data by these names (never by column position), so the
# order of columns in the user's data does not matter
BRGDGT_COLS <- c("fIa", "fIb", "fIc", "fIIa", "fIIa_", "fIIb", "fIIb_",
                 "fIIc", "fIIc_", "fIIIa", "fIIIa_", "fIIIb", "fIIIb_",
                 "fIIIc", "fIIIc_")

# Compounds grouped by the number of methyl groups on the alkyl chains:
# I = tetramethylated, II = pentamethylated, III = hexamethylated
# (a trailing underscore denotes the 6-methyl isomer, e.g. IIa')
METHYLATION_GROUPS <- list(
  tetra = c("fIa", "fIb", "fIc"),
  penta = c("fIIa", "fIIa_", "fIIb", "fIIb_", "fIIc", "fIIc_"),
  hexa = c("fIIIa", "fIIIa_", "fIIIb", "fIIIb_", "fIIIc", "fIIIc_")
)

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

#' Summarise brGDGTs by methylation group
#'
#' @description Sums the fractional abundances of tetra-, penta- and
#'   hexamethylated brGDGTs for each sample.
#'
#' @param data A data frame with the 15 brGDGT columns (fIa, fIb, fIc, fIIa,
#'   fIIa_, fIIb, fIIb_, fIIc, fIIc_, fIIIa, fIIIa_, fIIIb, fIIIb_, fIIIc,
#'   fIIIc_), as fractions or percentages.
#' @param id_col Name of the column holding sample identifiers. If absent,
#'   samples are numbered by row.
#'
#' @return A data frame with columns `id`, `tetra`, `penta` and `hexa`; each
#'   row sums to 1.
#' @export
#'
#' @examples
#' data(example_data)
#' head(brgdgt_methylation(example_data))
brgdgt_methylation <- function(data, id_col = "ID") {
  fractions <- brgdgt_fractions(data)
  data.frame(id = get_ids(data, id_col), methylation_groups(fractions))
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

# Training dataset, cached after the first load. An environment is used
# because, unlike a plain variable, it can be modified from inside a function
# after the package namespace has been locked.
.brgdgt_cache <- new.env(parent = emptyenv())

load_reference_data <- function() {
  # reading the rds takes a moment, so only do it once per session
  if (is.null(.brgdgt_cache$reference)) {
    .brgdgt_cache$reference <- readRDS(extdata_path("reference_data.rds"))
  }
  .brgdgt_cache$reference
}

# Path to a file shipped in inst/extdata. system.file() returns "" rather than
# an error when the file is missing, so turn that into an informative error.
extdata_path <- function(file) {
  path <- system.file("extdata", file, package = "envirolipids")
  if (!nzchar(path)) {
    stop("File '", file, "' not found in the envirolipids installation",
         call. = FALSE)
  }
  path
}

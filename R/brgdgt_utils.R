#' @importFrom ggplot2 autoplot .data
#' @importFrom randomForest randomForest
NULL

#' @export
ggplot2::autoplot

# The 15 brGDGT compounds used by the models
BRGDGT_COLS <- c("fIa", "fIb", "fIc", "fIIa", "fIIa_", "fIIb", "fIIb_",
                 "fIIc", "fIIc_", "fIIIa", "fIIIa_", "fIIIb", "fIIIb_",
                 "fIIIc", "fIIIc_")

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
  missing_cols <- setdiff(BRGDGT_COLS, names(data))
  if (length(missing_cols) > 0) {
    stop("`data` is missing brGDGT columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
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
  log_comp <- log(fractions + 1e-6)
  as.data.frame(log_comp - rowMeans(log_comp))
}

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

methylation_groups <- function(fractions) {
  groups <- vapply(METHYLATION_GROUPS,
                   function(cols) rowSums(fractions[, cols, drop = FALSE]),
                   numeric(nrow(fractions)))
  groups <- matrix(groups, nrow = nrow(fractions),
                   dimnames = list(NULL, names(METHYLATION_GROUPS)))
  as.data.frame(groups / rowSums(groups))
}

# Training dataset, cached after the first load
.brgdgt_cache <- new.env(parent = emptyenv())

load_reference_data <- function() {
  if (is.null(.brgdgt_cache$reference)) {
    .brgdgt_cache$reference <- readRDS(extdata_path("reference_data.rds"))
  }
  .brgdgt_cache$reference
}

extdata_path <- function(file) {
  path <- system.file("extdata", file, package = "envirolipids")
  if (!nzchar(path)) {
    stop("File '", file, "' not found in the envirolipids installation",
         call. = FALSE)
  }
  path
}

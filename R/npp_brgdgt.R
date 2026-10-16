#' Predict NPP from brGDGT data
#'
#' @description Predicts Net Primary Productivity (NPP) from the fractional
#'   abundances of 15 brGDGTs, using a centred log-ratio (CLR) transformation
#'   and a pre-trained random forest model. Check first whether your samples
#'   fall within the training data with [npp_brgdgt_applicable()].
#'
#' @param data A data frame with the 15 brGDGT columns (fIa, fIb, fIc, fIIa,
#'   fIIa_, fIIb, fIIb_, fIIc, fIIc_, fIIIa, fIIIa_, fIIIb, fIIIb_, fIIIc,
#'   fIIIc_), given either as fractions (summing to 1) or percentages
#'   (summing to 100). Columns are matched by name.
#' @param sample_type Sample type, either a single value or one per sample.
#'   One of "Lacustrine Sediment", "Lacustrine SPM", "Low DO Lacustrine SPM",
#'   "Peat", "Riverine Sediment and SPM", "Soil" or
#'   "Lacustrine Meso/Microcosm".
#' @param ci Logical; if `TRUE` (default), adds the 2.5% and 97.5% quantiles
#'   of the predictions across trees.
#' @param id_col Name of the column holding sample identifiers. If absent,
#'   samples are numbered by row.
#' @param model The fitted model; defaults to [load_npp_model()].
#'
#' @return A `brgdgt_prediction` data frame with columns `id`, `sample_type`,
#'   `estimate` and, if `ci = TRUE`, `lower` and `upper`. NPP is in
#'   gC m-2 yr-1. Use [ggplot2::autoplot()] to plot it and [hanpp()] to compare
#'   it with climate-based NPP.
#' @export
#'
#' @examples
#' \dontrun{
#' data(example_data)
#' npp <- npp_brgdgt(example_data, sample_type = "Lacustrine Sediment")
#' npp
#' ggplot2::autoplot(npp, x = example_data$Age)
#' }
npp_brgdgt <- function(data, sample_type, ci = TRUE, id_col = "ID",
                       model = load_npp_model()) {
  predict_rf_brgdgt(data, sample_type, model = model, ci = ci,
                    id_col = id_col, response = "NPP",
                    units = "gC m-2 yr-1")
}

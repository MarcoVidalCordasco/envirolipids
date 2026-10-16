#' Predict temperature from brGDGT data
#'
#' @description Predicts mean annual air temperature (MAT) from the fractional
#'   abundances of 15 brGDGTs, using a centred log-ratio (CLR) transformation
#'   and a pre-trained random forest model.
#'
#'   The temperature model is not yet shipped with the package; until it is,
#'   a fitted model must be supplied via `model`.
#'
#' @inheritParams npp_brgdgt
#' @param model The fitted model; defaults to [load_temp_model()].
#'
#' @return A `brgdgt_prediction` data frame with columns `id`, `sample_type`,
#'   `estimate` and, if `ci = TRUE`, `lower` and `upper`. MAT is in degrees C.
#' @export
#'
#' @examples
#' \dontrun{
#' data(example_data)
#' mat <- temp_brgdgt(example_data, sample_type = "Lacustrine Sediment")
#' ggplot2::autoplot(mat, x = example_data$Age)
#' }
temp_brgdgt <- function(data, sample_type, ci = TRUE, id_col = "ID",
                        model = load_temp_model()) {
  predict_rf_brgdgt(data, sample_type, model = model, ci = ci,
                    id_col = id_col, response = "MAT",
                    units = "\u00b0C")
}

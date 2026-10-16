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
  # Same pipeline as npp_brgdgt(), with the temperature model. The default
  # `model` is only evaluated when it is first used, so until the temperature
  # model ships, supplying `model` avoids load_temp_model()'s error.
  # The units use the escape u00b0 for the degree sign (written with a
  # backslash below), keeping the R code ASCII as CRAN requires.
  predict_rf_brgdgt(data, sample_type, model = model, ci = ci,
                    id_col = id_col, response = "MAT",
                    units = "\u00b0C")
}

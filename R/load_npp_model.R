#' Load the pre-trained brGDGT models
#'
#' @description `load_npp_model()` loads the random forest used by
#'   [npp_brgdgt()], trained on a global compilation of brGDGT data.
#'   `load_temp_model()` loads the model used by [temp_brgdgt()].
#'
#' @return A `train` object (caret) wrapping a `randomForest` model.
#' @export
#'
#' @examples
#' \dontrun{
#' model <- load_npp_model()
#' model$finalModel
#' }
load_npp_model <- function() {
  readRDS(extdata_path("model_brGDGT-NPP_all1.rds"))
}

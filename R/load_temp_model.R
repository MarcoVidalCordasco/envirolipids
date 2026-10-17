#' @rdname load_npp_model
#' @export
load_temp_model <- function() {
  # The temperature model has not been added to the package yet. Check for it
  # here, rather than through extdata_path(), to give a message explaining
  # how to work around its absence.
  path <- system.file("extdata", "model_brGDGT-MAT.rds",
                      package = "envirolipids")
  if (!nzchar(path)) {
    stop("The brGDGT temperature model is not yet included in envirolipids. ",
         "Pass a fitted model to `temp_brgdgt()` via `model`.", call. = FALSE)
  }
  readRDS(path)
}

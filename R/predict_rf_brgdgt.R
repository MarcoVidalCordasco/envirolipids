# Shared prediction code for the brGDGT random forest models.
# The models take as predictors the CLR-transformed compounds plus one 0/1
# dummy variable per sample type, so we build exactly that table for the new
# samples and pass it to the random forest.
predict_rf_brgdgt <- function(data, sample_type, model, ci, id_col,
                              response, units) {
  # CLR-transformed compounds: the first 15 predictors
  newdata <- brgdgt_clr(brgdgt_fractions(data))
  n <- nrow(newdata)

  # sample_type is either one value for all samples, or one per sample
  if (!is.character(sample_type) || !length(sample_type) %in% c(1, n)) {
    stop("`sample_type` must be a character vector of length 1 or ",
         "nrow(data)", call. = FALSE)
  }
  unknown <- setdiff(sample_type, SAMPLE_TYPES)
  if (length(unknown) > 0) {
    stop("Unknown `sample_type`: ", paste(unknown, collapse = ", "),
         ". Must be one of: ", paste(SAMPLE_TYPES, collapse = ", "),
         call. = FALSE)
  }
  sample_type <- rep_len(sample_type, n)
  # pool the aquatic SPM types, exactly as was done before training
  type_fixed <- ifelse(sample_type %in% AQUATIC_SPM_TYPES,
                       "Aquatic_SPM_combined", sample_type)

  # The saved model is a caret `train` object; the random forest itself is in
  # $finalModel. A bare randomForest can also be passed via `model`.
  rf <- if (inherits(model, "train")) model$finalModel else model
  # The model's own predictor names tell us which dummies it expects. For the
  # NPP model these are Sampletype_fixedLacustrine Sediment, ...Peat and
  # ...Soil: aquatic SPM is the baseline level, so it has no dummy of its own
  # and is coded as 0 for all three. Reading the names from the model, rather
  # than hard-coding them, means a model with different levels still works.
  predictors <- rownames(rf$importance)
  for (dummy in grep("^Sampletype_fixed", predictors, value = TRUE)) {
    # 1 for samples of this type, 0 otherwise
    newdata[[dummy]] <- as.numeric(type_fixed ==
                                     sub("^Sampletype_fixed", "", dummy))
  }
  # any predictor we could not build means the model was trained differently
  missing_predictors <- setdiff(predictors, names(newdata))
  if (length(missing_predictors) > 0) {
    stop("Model uses unexpected predictors: ",
         paste(missing_predictors, collapse = ", "), call. = FALSE)
  }
  # put the columns in the order the model was trained with
  newdata <- newdata[, predictors, drop = FALSE]

  out <- data.frame(id = get_ids(data, id_col), sample_type = sample_type)
  if (ci) {
    # predict.all = TRUE also returns each tree's prediction ($individual, a
    # samples x trees matrix) besides the forest mean ($aggregate). The spread
    # across trees gives a 95% interval. Prediction is deterministic, so the
    # interval is the same every time.
    preds <- stats::predict(rf, newdata = newdata, predict.all = TRUE)
    # apply() over rows returns a 2 x samples matrix: row 1 = 2.5%,
    # row 2 = 97.5%
    bounds <- apply(preds$individual, 1, stats::quantile,
                    probs = c(0.025, 0.975), na.rm = TRUE)
    out$estimate <- preds$aggregate
    out$lower <- bounds[1, ]
    out$upper <- bounds[2, ]
  } else {
    out$estimate <- stats::predict(rf, newdata = newdata)
  }
  rownames(out) <- NULL

  # A data frame subclass, so it can be used as a normal table. The
  # attributes record what was predicted, for printing and plotting.
  structure(out, class = c("brgdgt_prediction", "data.frame"),
            response = response, units = units)
}

# Print method for the brgdgt_prediction class returned by predict_rf_brgdgt()
#' @export
print.brgdgt_prediction <- function(x, ...) {
  cat("brGDGT-based ", attr(x, "response"), " predictions (",
      attr(x, "units"), ") for ", nrow(x), " samples\n", sep = "")
  if ("lower" %in% names(x)) {
    cat("lower/upper: 2.5% and 97.5% quantiles across random forest trees\n")
  }
  cat("\n")
  # print as a plain data frame, to avoid calling this method recursively
  print(as.data.frame(x), ...)
  invisible(x)
}

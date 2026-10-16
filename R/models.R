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

#' @rdname load_npp_model
#' @export
load_temp_model <- function() {
  path <- system.file("extdata", "model_brGDGT-MAT.rds", package = "envirolipids")
  if (!nzchar(path)) {
    stop("The brGDGT temperature model is not yet included in envirolipids. ",
         "Pass a fitted model to `temp_brgdgt()` via `model`.", call. = FALSE)
  }
  readRDS(path)
}

SAMPLE_TYPES <- c("Lacustrine Sediment", "Lacustrine SPM",
                  "Low DO Lacustrine SPM", "Peat",
                  "Riverine Sediment and SPM", "Soil",
                  "Lacustrine Meso/Microcosm")

# Sample types pooled into a single class when training the models
AQUATIC_SPM_TYPES <- c("Lacustrine SPM", "Low DO Lacustrine SPM",
                       "Riverine Sediment and SPM", "Lacustrine Meso/Microcosm")

# Shared prediction code for the brGDGT random forest models
predict_rf_brgdgt <- function(data, sample_type, model, ci, id_col,
                              response, units) {
  newdata <- brgdgt_clr(brgdgt_fractions(data))
  n <- nrow(newdata)

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
  type_fixed <- ifelse(sample_type %in% AQUATIC_SPM_TYPES,
                       "Aquatic_SPM_combined", sample_type)

  rf <- if (inherits(model, "train")) model$finalModel else model
  predictors <- rownames(rf$importance)
  for (dummy in grep("^Sampletype_fixed", predictors, value = TRUE)) {
    newdata[[dummy]] <- as.numeric(type_fixed ==
                                     sub("^Sampletype_fixed", "", dummy))
  }
  missing_predictors <- setdiff(predictors, names(newdata))
  if (length(missing_predictors) > 0) {
    stop("Model uses unexpected predictors: ",
         paste(missing_predictors, collapse = ", "), call. = FALSE)
  }
  newdata <- newdata[, predictors, drop = FALSE]

  out <- data.frame(id = get_ids(data, id_col), sample_type = sample_type)
  if (ci) {
    preds <- stats::predict(rf, newdata = newdata, predict.all = TRUE)
    bounds <- apply(preds$individual, 1, stats::quantile,
                    probs = c(0.025, 0.975), na.rm = TRUE)
    out$estimate <- preds$aggregate
    out$lower <- bounds[1, ]
    out$upper <- bounds[2, ]
  } else {
    out$estimate <- stats::predict(rf, newdata = newdata)
  }
  rownames(out) <- NULL

  structure(out, class = c("brgdgt_prediction", "data.frame"),
            response = response, units = units)
}

#' @export
print.brgdgt_prediction <- function(x, ...) {
  cat("brGDGT-based ", attr(x, "response"), " predictions (",
      attr(x, "units"), ") for ", nrow(x), " samples\n", sep = "")
  if ("lower" %in% names(x)) {
    cat("lower/upper: 2.5% and 97.5% quantiles across random forest trees\n")
  }
  cat("\n")
  print(as.data.frame(x), ...)
  invisible(x)
}

#' Plot brGDGT-based predictions
#'
#' @description Plots predicted values (with the 95% interval across trees
#'   if available) against sample id or a supplied variable such as age.
#'
#' @param object A `brgdgt_prediction` object from [npp_brgdgt()] or
#'   [temp_brgdgt()].
#' @param x Optional numeric vector (e.g. sample age or depth), one value per
#'   sample, to use as the x axis. Defaults to the sample id.
#' @param xlab Label for the x axis.
#' @param ... Not used.
#'
#' @return A ggplot object.
#' @export
autoplot.brgdgt_prediction <- function(object, x = NULL, xlab = NULL, ...) {
  df <- as.data.frame(object)
  df$x <- plot_x(df, x)
  if (is.null(xlab)) xlab <- if (is.null(x)) "Sample" else deparse(substitute(x))
  ylab <- paste0(attr(object, "response"), " (", attr(object, "units"), ")")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$estimate))
  if ("lower" %in% names(df)) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      fill = "grey80", alpha = 0.6)
  }
  p + ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(x = xlab, y = ylab,
                  title = paste("brGDGT-based", attr(object, "response"))) +
    ggplot2::theme_minimal()
}

# x axis for series plots: a user-supplied vector, else the sample id
plot_x <- function(df, x) {
  if (is.null(x)) {
    return(if (is.numeric(df$id)) df$id else seq_len(nrow(df)))
  }
  if (length(x) != nrow(df)) {
    stop("`x` must have one value per sample", call. = FALSE)
  }
  x
}

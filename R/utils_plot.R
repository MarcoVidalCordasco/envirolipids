# Internal helpers shared by the autoplot() methods.

# x values for the plots: a factor of sample labels (plotted as separate
# samples) or a numeric vector (plotted as a continuous series)
plot_x <- function(df, x) {
  if (is.null(x)) {
    # By default, plot by sample id. make.unique() guards against repeated
    # ids, which would otherwise be merged onto one axis position. The factor
    # levels keep the samples in the order of the input data.
    labels <- make.unique(as.character(df$id))
    return(factor(labels, levels = labels))
  }
  if (length(x) != nrow(df)) {
    stop("`x` must have one value per sample", call. = FALSE)
  }
  if (is.numeric(x)) x else factor(x, levels = unique(x))
}

# For plots by sample, angle the sample labels so that long ids do not overlap
theme_sample_axis <- function(x) {
  if (is.factor(x)) {
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
}

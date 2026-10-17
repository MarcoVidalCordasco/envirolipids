#' Plot brGDGT-based predictions
#'
#' @description Plots predicted values, with the 95% interval across trees if
#'   available. There are two layouts:
#'   * by sample (the default, e.g. for modern samples): one point per
#'     sample, with error bars, labelled by sample id;
#'   * along a numeric `x` (e.g. age or depth for a sediment core): a line
#'     through the samples, with the interval as a ribbon.
#'
#' @param object A `brgdgt_prediction` object from [npp_brgdgt()] or
#'   [temp_brgdgt()].
#' @param x Optional vector with one value per sample. A numeric vector (e.g.
#'   sample age or depth) gives a continuous axis; anything else is treated as
#'   sample labels. If `NULL` (default), samples are plotted by id.
#' @param xlab Label for the x axis.
#' @param ... Not used.
#'
#' @return A ggplot object.
#' @export
autoplot.brgdgt_prediction <- function(object, x = NULL, xlab = NULL, ...) {
  df <- as.data.frame(object)
  df$x <- plot_x(df, x)
  # default x label: the expression the user passed (e.g. "example_data$Age")
  if (is.null(xlab)) {
    xlab <- if (is.null(x)) "Sample" else deparse(substitute(x))
  }
  ylab <- paste0(attr(object, "response"), " (", attr(object, "units"), ")")
  # the interval only exists if the predictions were made with ci = TRUE
  has_ci <- "lower" %in% names(df)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$estimate))
  if (is.factor(df$x)) {
    # separate samples: no line between them, and error bars for the interval
    if (has_ci) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        width = 0.3, colour = "grey50"
      )
    }
    p <- p + ggplot2::geom_point(size = 2)
  } else {
    # a series (e.g. through time): connect the samples, interval as a ribbon
    if (has_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        fill = "grey80", alpha = 0.6
      )
    }
    p <- p + ggplot2::geom_line() + ggplot2::geom_point()
  }
  p + ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_minimal() +
    theme_sample_axis(df$x)
}

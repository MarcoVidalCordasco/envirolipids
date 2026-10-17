#' Plot climate-based NPP
#'
#' @param object An `npp_climate` object from [npp_miami()] or
#'   [npp_thornthwaite()].
#' @inheritParams autoplot.brgdgt_prediction
#'
#' @return A ggplot object. For the Miami model, points are coloured by the
#'   limiting factor.
#' @export
autoplot.npp_climate <- function(object, x = NULL, xlab = NULL, ...) {
  df <- as.data.frame(object)
  df$x <- plot_x(df, x)
  if (is.null(xlab)) {
    xlab <- if (is.null(x)) "Sample" else deparse(substitute(x))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$npp))
  # connect the points only for a continuous series (e.g. through time), not
  # between separate samples
  if (!is.factor(df$x)) {
    p <- p + ggplot2::geom_line(colour = "grey50")
  }
  # only the Miami model reports a limiting factor, so colour by it if present
  if ("limiting" %in% names(df)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$limiting)) +
      ggplot2::labs(colour = "Limiting factor")
  } else {
    p <- p + ggplot2::geom_point()
  }
  p + ggplot2::labs(x = xlab, y = paste0("NPP (", attr(object, "units"), ")")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    theme_sample_axis(df$x)
}

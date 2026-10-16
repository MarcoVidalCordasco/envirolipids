#' Compute the deficit between climate-based and empirical NPP
#'
#' @description Computes the difference between potential NPP expected from
#'   climate and empirical (brGDGT-based) NPP. Positive values mean NPP is
#'   lower than climate predicts, as expected under human appropriation of
#'   NPP (HANPP); negative values mean NPP is higher than expected.
#'
#' @param npp_empirical Empirical NPP: a `brgdgt_prediction` from
#'   [npp_brgdgt()] or a numeric vector.
#' @param npp_model Climate-based NPP: an `npp_climate` object from
#'   [npp_miami()] or [npp_thornthwaite()], or a numeric vector of the same
#'   length as `npp_empirical`.
#'
#' @return A `hanpp` data frame with columns `id`, `npp_empirical`,
#'   `npp_model` and `hanpp` (`npp_model - npp_empirical`). If
#'   `npp_empirical` has an interval (`lower`, `upper`), `hanpp_lower` and
#'   `hanpp_upper` are added.
#' @export
#'
#' @examples
#' hanpp(npp_empirical = c(550, 680, 720),
#'       npp_model = npp_miami(mat = c(8, 10, 12), map = c(600, 800, 900)))
hanpp <- function(npp_empirical, npp_model) {
  # reduce both inputs, whatever their type, to the same simple list of
  # values (plus ids and interval when available)
  emp <- extract_npp(npp_empirical, "npp_empirical")
  mod <- extract_npp(npp_model, "npp_model")
  if (length(emp$npp) != length(mod$npp)) {
    stop("`npp_empirical` and `npp_model` must have the same number of ",
         "samples", call. = FALSE)
  }
  # Ids come from the brGDGT predictions if available, then from the climate
  # model, and otherwise default to sample numbers. The two inputs are matched
  # by position, not by id, so they must be in the same sample order.
  id <- if (!is.null(emp$id)) emp$id else if (!is.null(mod$id)) mod$id else
    seq_along(emp$npp)

  out <- data.frame(id = id, npp_empirical = emp$npp, npp_model = mod$npp,
                    hanpp = mod$npp - emp$npp)
  # Carry the brGDGT prediction interval through to the deficit. Because the
  # deficit is model - empirical, the bounds swap: the highest empirical NPP
  # gives the smallest deficit. Uncertainty in the climate model is not
  # included, as the climate models give no interval.
  if (!is.null(emp$lower)) {
    out$hanpp_lower <- mod$npp - emp$upper
    out$hanpp_upper <- mod$npp - emp$lower
  }
  structure(out, class = c("hanpp", "data.frame"), model = mod$model)
}

# Turn any accepted input into list(npp, id, lower, upper, model); elements
# that do not apply are left NULL. `arg` is the argument name, used in error
# messages.
extract_npp <- function(x, arg) {
  if (inherits(x, "brgdgt_prediction")) {
    # guard against passing temp_brgdgt() output by mistake: it has the same
    # class but holds temperatures
    if (attr(x, "response") != "NPP") {
      stop("`", arg, "` holds ", attr(x, "response"), " predictions, not NPP",
           call. = FALSE)
    }
    return(list(npp = x$estimate, id = x$id, lower = x$lower, upper = x$upper))
  }
  if (inherits(x, "npp_climate")) {
    return(list(npp = x$npp, id = x$id, model = attr(x, "model")))
  }
  # as.vector() drops names and dimensions, e.g. from a one-column matrix
  if (is.numeric(x)) {
    return(list(npp = as.vector(x)))
  }
  stop("`", arg, "` must be a numeric vector or the output of npp_brgdgt(), ",
       "npp_miami() or npp_thornthwaite()", call. = FALSE)
}

#' @export
print.hanpp <- function(x, ...) {
  model <- attr(x, "model")
  cat("NPP deficit (climate", if (!is.null(model)) paste0("[", model, "] "),
      "- empirical, gC m-2 yr-1) for ", nrow(x), " samples\n",
      "Mean deficit: ", round(mean(x$hanpp, na.rm = TRUE), 1), "\n\n",
      sep = "")
  print(as.data.frame(x), ...)
  invisible(x)
}

#' Plot empirical and climate-based NPP and their difference
#'
#' @param object A `hanpp` object from [hanpp()].
#' @inheritParams autoplot.brgdgt_prediction
#'
#' @return A ggplot object showing empirical NPP, climate-based NPP and the
#'   deficit, with the deficit interval if available.
#' @export
autoplot.hanpp <- function(object, x = NULL, xlab = NULL, ...) {
  df <- as.data.frame(object)
  df$x <- plot_x(df, x)
  if (is.null(xlab)) xlab <- if (is.null(x)) "Sample" else deparse(substitute(x))

  # Reshape to long format (one row per sample and series), so that ggplot
  # draws the three series as coloured lines with a single legend. `series`
  # maps column names to legend labels.
  series <- c(npp_empirical = "Empirical NPP", npp_model = "Climate NPP",
              hanpp = "Deficit (HANPP)")
  long <- do.call(rbind, lapply(names(series), function(col) {
    data.frame(x = df$x, value = df[[col]], series = series[[col]])
  }))
  # fix the legend order to the order above, not alphabetical
  long$series <- factor(long$series, levels = series)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$x, y = .data$value,
                                          colour = .data$series))
  # The interval band around the deficit uses the wide data frame, so it must
  # not inherit the long-format colour mapping
  if ("hanpp_lower" %in% names(df)) {
    p <- p + ggplot2::geom_ribbon(
      data = df,
      ggplot2::aes(x = .data$x, ymin = .data$hanpp_lower,
                   ymax = .data$hanpp_upper),
      inherit.aes = FALSE, fill = "grey80", alpha = 0.6)
  }
  # zero line: deficits above it mean NPP is lower than climate predicts
  p + ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_manual(values = c("#D55E00", "#0072B2", "black")) +
    ggplot2::labs(x = xlab, y = "NPP (gC m-2 yr-1)", colour = NULL,
                  title = "Empirical vs climate-based NPP") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

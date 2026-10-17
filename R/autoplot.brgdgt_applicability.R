#' Diagnostic plots for the applicability check
#'
#' @description Plots the results of [npp_brgdgt_applicable()]:
#'   * `"mahalanobis"`: Mahalanobis distance of each sample, with the cut-off,
#'     by sample or along a numeric `x` such as age;
#'   * `"pca"`: samples on the first two principal components, over the
#'     training data;
#'   * `"ternary"`: ternary plot of methylation groups, with training samples
#'     sized by NPP;
#'   * `"map"`: locations of the training sites most similar in composition
#'     to the samples (requires the maps package).
#'
#' @param object A `brgdgt_applicability` object.
#' @param type Which plot to draw.
#' @param x For `type = "mahalanobis"`, an optional vector with one value per
#'   sample to use as the x axis (e.g. age or depth). If `NULL` (default),
#'   samples are plotted by id.
#' @param xlab For `type = "mahalanobis"`, label for the x axis.
#' @param n_sites For `type = "map"`, number of most similar sites to show.
#' @param ... Not used.
#'
#' @return A ggplot object.
#' @export
autoplot.brgdgt_applicability <- function(object,
                                          type = c("mahalanobis", "pca",
                                                   "ternary", "map"),
                                          x = NULL, xlab = NULL,
                                          n_sites = 10, ...) {
  # match.arg() picks the first type by default and gives a helpful error for
  # unknown types. Each type has its own plotting function below, which only
  # reads from `object` and never recomputes the diagnostics.
  type <- match.arg(type)
  if (is.null(xlab)) {
    xlab <- if (is.null(x)) "Sample" else deparse(substitute(x))
  }
  switch(
    type,
    mahalanobis = plot_mahalanobis(object, x, xlab = xlab),
    pca = plot_pca(object),
    ternary = plot_ternary(object),
    map = plot_similarity_map(object, n_sites)
  )
}

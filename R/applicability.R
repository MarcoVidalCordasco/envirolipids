#' Check whether brGDGT samples fall within the model's training data
#'
#' @description Assesses whether new brGDGT samples fall within the
#'   composition space of the data used to train the NPP model, which helps
#'   flag samples whose predictions may be unreliable. Compositions are
#'   CLR-transformed, and for each sample this computes:
#'   * the Mahalanobis distance to the training data (on the principal
#'     components of the training data explaining 95% of its variance if
#'     `use_pca = TRUE`, otherwise on a regularised CLR covariance), with a
#'     chi-squared p-value. New samples are projected onto the training PCA,
#'     so each sample's result does not depend on the others in `data`;
#'   * the Euclidean distances to its nearest training samples.
#'
#'   Use `print()` for a table of results and [ggplot2::autoplot()] for
#'   diagnostic plots.
#'
#' @inheritParams npp_brgdgt
#' @param threshold Quantile of the chi-squared distribution used as the
#'   Mahalanobis distance cut-off for being within the domain (default 0.95).
#' @param use_pca Logical; if `TRUE` (default), computes Mahalanobis distances
#'   on principal components to avoid a singular covariance matrix.
#' @param n_neighbors Number of nearest training samples to find for each
#'   sample (default 5).
#'
#' @return A `brgdgt_applicability` object, a list with elements:
#'   \item{samples}{Data frame with, for each sample, `id`, `mahalanobis`,
#'     `p_value`, `within_domain` and `nn_distance` (distance to the nearest
#'     training sample).}
#'   \item{neighbours}{Data frame of nearest training samples for each sample:
#'     `id`, `rank`, `ref_index`, `distance`, `site`, `longitude`,
#'     `latitude`, `npp`.}
#'   \item{threshold_value}{Mahalanobis distance cut-off.}
#'   \item{method}{`"pca"` or `"regularized"`.}
#'   \item{n_pcs}{Degrees of freedom (number of PCs or CLR components).}
#'   \item{pca}{PC1 and PC2 scores of the training (`ref`) and new (`new`)
#'     samples, and the proportion of variance each explains (`var_explained`).}
#'   \item{methylation}{Methylation group proportions of the training (`ref`)
#'     and new (`new`) samples, see [brgdgt_methylation()].}
#'   \item{reference}{Training data NPP, sites and coordinates.}
#' @export
#'
#' @examples
#' \dontrun{
#' data(example_data)
#' check <- npp_brgdgt_applicable(example_data)
#' check
#' ggplot2::autoplot(check, type = "pca")
#' }
npp_brgdgt_applicable <- function(data, threshold = 0.95, use_pca = TRUE,
                                  n_neighbors = 5, id_col = "ID") {
  reference <- load_reference_data()
  ids <- get_ids(data, id_col)
  fractions <- brgdgt_fractions(data)
  clr_new <- as.matrix(brgdgt_clr(fractions))
  clr_ref <- as.matrix(reference$clr_values[, BRGDGT_COLS])
  n_new <- nrow(clr_new)

  # Euclidean distances in CLR space, new samples x training samples
  dists <- sqrt(pmax(outer(rowSums(clr_new^2), rowSums(clr_ref^2), "+") -
                       2 * tcrossprod(clr_new, clr_ref), 0))
  n_neighbors <- min(n_neighbors, ncol(dists))
  neighbours <- do.call(rbind, lapply(seq_len(n_new), function(i) {
    nn <- order(dists[i, ])[seq_len(n_neighbors)]
    data.frame(id = ids[i], rank = seq_len(n_neighbors), ref_index = nn,
               distance = dists[i, nn], site = reference$sites[nn],
               longitude = reference$coordinates$Longitude[nn],
               latitude = reference$coordinates$Latitude[nn],
               npp = reference$npp[nn])
  }))

  # PCA fitted on the training data, with new samples projected onto it
  pca <- stats::prcomp(clr_ref, center = TRUE, scale. = TRUE)
  var_explained <- pca$sdev^2 / sum(pca$sdev^2)
  scores_ref <- pca$x
  scores_new <- stats::predict(pca, newdata = clr_new)

  if (use_pca) {
    n_pcs <- min(which(cumsum(var_explained) >= 0.95))
    space_ref <- scores_ref[, seq_len(n_pcs), drop = FALSE]
    space_new <- scores_new[, seq_len(n_pcs), drop = FALSE]
    method <- "pca"
  } else {
    n_pcs <- ncol(clr_ref)
    space_ref <- clr_ref
    space_new <- clr_new
    method <- "regularized"
  }
  cov_matrix <- stats::cov(space_ref) + diag(1e-6, n_pcs)
  mahal <- stats::mahalanobis(space_new, colMeans(space_ref), cov_matrix)
  threshold_value <- stats::qchisq(threshold, df = n_pcs)

  samples <- data.frame(
    id = ids,
    mahalanobis = unname(mahal),
    p_value = stats::pchisq(mahal, df = n_pcs, lower.tail = FALSE),
    within_domain = unname(mahal <= threshold_value),
    nn_distance = apply(dists, 1, min)
  )

  structure(
    list(
      samples = samples,
      neighbours = neighbours,
      threshold = threshold,
      threshold_value = threshold_value,
      method = method,
      n_pcs = n_pcs,
      pca = list(ref = scores_ref[, 1:2, drop = FALSE],
                 new = scores_new[, 1:2, drop = FALSE],
                 var_explained = var_explained[1:2]),
      methylation = list(ref = methylation_groups(reference$fractions),
                         new = methylation_groups(fractions)),
      reference = list(npp = reference$npp, sites = reference$sites,
                       coordinates = reference$coordinates)
    ),
    class = "brgdgt_applicability"
  )
}

#' @export
print.brgdgt_applicability <- function(x, ...) {
  s <- x$samples
  cat("brGDGT applicability domain check\n",
      "Mahalanobis distance (", x$method, ", df = ", x$n_pcs,
      "), cut-off ", round(x$threshold_value, 2), " (", x$threshold,
      " quantile)\n",
      "Samples within domain: ", sum(s$within_domain), " of ", nrow(s),
      " (", round(mean(s$within_domain) * 100, 1), "%)\n\n", sep = "")
  out <- s
  out$mahalanobis <- round(out$mahalanobis, 2)
  out$p_value <- format.pval(out$p_value, digits = 3)
  out$nn_distance <- round(out$nn_distance, 3)
  print(out, ...)
  invisible(x)
}

#' @export
as.data.frame.brgdgt_applicability <- function(x, ...) {
  x$samples
}

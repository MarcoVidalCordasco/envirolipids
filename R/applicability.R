#' Check whether brGDGT samples fall within the model's training data
#'
#' @description Assesses whether new brGDGT samples fall within the
#'   composition space of the data used to train the NPP model, which helps
#'   flag samples whose predictions may be unreliable. Compositions are
#'   CLR-transformed, and for each sample this computes:
#'   * the Mahalanobis distance to the training data (on the principal
#'     components of the training data explaining 95% of its variance if
#'     `use_pca = TRUE`, otherwise on a regularised CLR covariance). New
#'     samples are projected onto the training PCA, so each sample's result
#'     does not depend on the others in `data`;
#'   * the Euclidean distances to its nearest training samples.
#'
#'   Use `print()` for a table of results and [ggplot2::autoplot()] for
#'   diagnostic plots.
#'
#'   A sample is within the domain if its Mahalanobis distance is no larger
#'   than that of a given proportion (`threshold`) of the training samples.
#'   The cut-off is taken from the training samples' own distances rather
#'   than from a chi-squared distribution: the chi-squared assumes
#'   multivariate normal data, but brGDGT compositions are heavy-tailed, so it
#'   would flag a large share of the training data itself as out of domain.
#'   With the default `threshold = 0.99`, only the most extreme 1% of training
#'   samples lie beyond the cut-off.
#'
#' @inheritParams npp_brgdgt
#' @param threshold Proportion of the training samples that fall within the
#'   domain: the Mahalanobis distance cut-off is this quantile of the training
#'   samples' distances (default 0.99).
#' @param use_pca Logical; if `TRUE` (default), computes Mahalanobis distances
#'   on principal components to avoid a singular covariance matrix.
#' @param n_neighbors Number of nearest training samples to find for each
#'   sample (default 5).
#'
#' @return A `brgdgt_applicability` object, a list with elements:
#'   \item{samples}{Data frame with, for each sample, `id`, `mahalanobis`,
#'     `p_value` (the proportion of training samples at least as far from the
#'     centre of the training data), `within_domain` and `nn_distance`
#'     (distance to the nearest training sample).}
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
npp_brgdgt_applicable <- function(data, threshold = 0.99, use_pca = TRUE,
                                  n_neighbors = 5, id_col = "ID") {
  reference <- load_reference_data()
  ids <- get_ids(data, id_col)
  # Put new and training samples in the same CLR space as the models. They are
  # held as matrices, since the distance computations are matrix algebra.
  fractions <- brgdgt_fractions(data)
  clr_new <- as.matrix(brgdgt_clr(fractions))
  clr_ref <- as.matrix(reference$clr_values[, BRGDGT_COLS])
  n_new <- nrow(clr_new)

  # --- Nearest training samples ---------------------------------------------
  # Euclidean distances in CLR space, new samples x training samples. Rather
  # than looping over pairs, use |a - b|^2 = |a|^2 + |b|^2 - 2 a.b for all
  # pairs at once: outer() adds every |a|^2 to every |b|^2, and tcrossprod()
  # gives every dot product a.b. Rounding can push values that should be 0
  # slightly negative, so pmax() clips them before the square root.
  dists <- sqrt(pmax(outer(rowSums(clr_new^2), rowSums(clr_ref^2), "+") -
                       2 * tcrossprod(clr_new, clr_ref), 0))
  # cannot ask for more neighbours than there are training samples
  n_neighbors <- min(n_neighbors, ncol(dists))
  # For each new sample, the n closest training samples (nearest first), with
  # where they come from and their NPP. This feeds the similarity map and lets
  # users inspect what their samples resemble.
  neighbours <- do.call(rbind, lapply(seq_len(n_new), function(i) {
    nn <- order(dists[i, ])[seq_len(n_neighbors)]
    data.frame(id = ids[i], rank = seq_len(n_neighbors), ref_index = nn,
               distance = dists[i, nn], site = reference$sites[nn],
               longitude = reference$coordinates$Longitude[nn],
               latitude = reference$coordinates$Latitude[nn],
               npp = reference$npp[nn])
  }))

  # --- Mahalanobis distance to the training data ----------------------------
  # Mahalanobis distance measures how far a sample is from the centre of the
  # training data, allowing for the spread and correlation of the compounds.
  # The CLR covariance matrix is singular (each CLR row sums to zero), so it
  # cannot be inverted directly. By default we work instead on the principal
  # components that carry 95% of the variance.
  #
  # The PCA is fitted on the training data only, and new samples are projected
  # onto it with predict(), so the reference space is fixed: a sample's result
  # does not depend on which other samples are checked alongside it.
  pca <- stats::prcomp(clr_ref, center = TRUE, scale. = TRUE)
  var_explained <- pca$sdev^2 / sum(pca$sdev^2)
  scores_ref <- pca$x
  scores_new <- stats::predict(pca, newdata = clr_new)

  if (use_pca) {
    # smallest number of components explaining at least 95% of the variance
    n_pcs <- min(which(cumsum(var_explained) >= 0.95))
    space_ref <- scores_ref[, seq_len(n_pcs), drop = FALSE]
    space_new <- scores_new[, seq_len(n_pcs), drop = FALSE]
    method <- "pca"
  } else {
    # use the full CLR space; the ridge added below makes the covariance
    # invertible
    n_pcs <- ncol(clr_ref)
    space_ref <- clr_ref
    space_new <- clr_new
    method <- "regularized"
  }
  # The small ridge on the diagonal guards against (near-)singular covariance
  # matrices. In PCA space it has negligible effect; in CLR space it is what
  # makes the matrix invertible.
  cov_matrix <- stats::cov(space_ref) + diag(1e-6, n_pcs)
  # stats::mahalanobis() returns the *squared* distance.
  mahal <- stats::mahalanobis(space_new, colMeans(space_ref), cov_matrix)
  # The cut-off and p-values come from the distances of the training samples
  # themselves. For multivariate normal data the distances would follow a
  # chi-squared distribution, but CLR compositions are heavy-tailed: about 14%
  # of the training samples lie beyond the chi-squared 95% quantile, so that
  # cut-off rejects much of the data the model was trained on. The empirical
  # distribution makes no assumption about the shape of the data, and defines
  # the domain directly in terms of the training samples.
  mahal_ref <- stats::mahalanobis(space_ref, colMeans(space_ref), cov_matrix)
  threshold_value <- unname(stats::quantile(mahal_ref, threshold))

  samples <- data.frame(
    id = ids,
    mahalanobis = unname(mahal),
    # proportion of training samples at least this far out; small values
    # flag unusual compositions
    p_value = vapply(mahal, function(m) mean(mahal_ref >= m), numeric(1)),
    within_domain = unname(mahal <= threshold_value),
    nn_distance = apply(dists, 1, min)
  )

  # Everything the print and autoplot methods need is computed here and
  # stored in the object, so that the methods only draw and never recompute.
  # Only PC1 and PC2 are kept, since they are all the PCA plot shows.
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
      "), cut-off ", round(x$threshold_value, 2), " (", x$threshold * 100,
      "% of training samples within)\n",
      "Samples within domain: ", sum(s$within_domain), " of ", nrow(s),
      " (", round(mean(s$within_domain) * 100, 1), "%)\n\n", sep = "")
  # round for display only; the stored values keep full precision
  out <- s
  out$mahalanobis <- round(out$mahalanobis, 2)
  out$p_value <- format.pval(out$p_value, digits = 3)
  out$nn_distance <- round(out$nn_distance, 3)
  print(out, ...)
  invisible(x)
}

# The per-sample table is what most users will want to keep or export
#' @export
as.data.frame.brgdgt_applicability <- function(x, ...) {
  x$samples
}

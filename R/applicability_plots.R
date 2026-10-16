#' Diagnostic plots for the applicability check
#'
#' @description Plots the results of [npp_brgdgt_applicable()]:
#'   * `"mahalanobis"`: Mahalanobis distance of each sample, with the cut-off;
#'   * `"pca"`: samples on the first two principal components, over the
#'     training data;
#'   * `"ternary"`: ternary plot of methylation groups, with training samples
#'     sized by NPP;
#'   * `"map"`: locations of the training sites most similar in composition
#'     to the samples (requires the maps package).
#'
#' @param object A `brgdgt_applicability` object.
#' @param type Which plot to draw.
#' @param n_sites For `type = "map"`, number of most similar sites to show.
#' @param ... Not used.
#'
#' @return A ggplot object.
#' @export
autoplot.brgdgt_applicability <- function(object,
                                          type = c("mahalanobis", "pca",
                                                   "ternary", "map"),
                                          n_sites = 10, ...) {
  type <- match.arg(type)
  switch(type,
         mahalanobis = plot_mahalanobis(object),
         pca = plot_pca(object),
         ternary = plot_ternary(object),
         map = plot_similarity_map(object, n_sites))
}

DOMAIN_COLOURS <- c("TRUE" = "#2E8B57", "FALSE" = "#DC143C")

scale_domain <- function() {
  ggplot2::scale_colour_manual(values = DOMAIN_COLOURS,
                               labels = c("TRUE" = "Within",
                                          "FALSE" = "Outside"),
                               name = "Domain", drop = FALSE)
}

domain_factor <- function(x) factor(x, levels = c("TRUE", "FALSE"))

plot_mahalanobis <- function(object) {
  df <- object$samples
  df$index <- seq_len(nrow(df))
  df$within_domain <- domain_factor(df$within_domain)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$index, y = .data$mahalanobis,
                                   colour = .data$within_domain)) +
    ggplot2::geom_hline(yintercept = object$threshold_value,
                        linetype = "dashed", colour = "red", alpha = 0.7) +
    ggplot2::geom_point(size = 3) +
    scale_domain() +
    ggplot2::labs(title = "Mahalanobis distance to training data",
                  x = "Sample index", y = "Mahalanobis distance") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

plot_pca <- function(object) {
  ref <- as.data.frame(object$pca$ref)
  new <- as.data.frame(object$pca$new)
  new$within_domain <- domain_factor(object$samples$within_domain)
  pct <- round(object$pca$var_explained * 100, 1)
  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$PC1, y = .data$PC2)) +
    ggplot2::geom_point(data = ref, colour = "grey70", alpha = 0.4,
                        size = 1.5) +
    ggplot2::geom_point(data = new,
                        ggplot2::aes(colour = .data$within_domain),
                        size = 3) +
    scale_domain() +
    ggplot2::labs(title = "Samples vs training data in CLR space",
                  x = paste0("PC1 (", pct[1], "%)"),
                  y = paste0("PC2 (", pct[2], "%)")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

# Ternary coordinates: tetra at the top, penta bottom left, hexa bottom right
ternary_xy <- function(df) {
  data.frame(x = df$hexa + df$tetra / 2, y = df$tetra * sqrt(3) / 2)
}

plot_ternary <- function(object) {
  ref <- ternary_xy(object$methylation$ref)
  ref$npp <- object$reference$npp
  new <- ternary_xy(object$methylation$new)
  new$within_domain <- domain_factor(object$samples$within_domain)

  h <- sqrt(3) / 2
  frame <- data.frame(x = c(0, 1, 0.5, 0), y = c(0, 0, h, 0))
  # grid lines at 20% steps, parallel to each side
  g <- seq(0.2, 0.8, by = 0.2)
  grid <- rbind(
    data.frame(x = g / 2, y = g * h, xend = 1 - g / 2, yend = g * h),
    data.frame(x = g, y = 0, xend = g + (1 - g) / 2, yend = (1 - g) * h),
    data.frame(x = 1 - g, y = 0, xend = (1 - g) / 2, yend = (1 - g) * h)
  )
  labels <- data.frame(x = c(0.5, -0.03, 1.03), y = c(h + 0.05, -0.04, -0.04),
                       label = c("Tetra", "Penta", "Hexa"),
                       hjust = c(0.5, 1, 0))

  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_segment(data = grid,
                          ggplot2::aes(xend = .data$xend, yend = .data$yend),
                          colour = "grey90") +
    ggplot2::geom_path(data = frame, colour = "grey40") +
    # training samples without NPP are drawn at a fixed size
    ggplot2::geom_point(data = ref[is.na(ref$npp), ], colour = "grey60",
                        alpha = 0.3, size = 0.3) +
    ggplot2::geom_point(data = ref[!is.na(ref$npp), ],
                        ggplot2::aes(size = .data$npp),
                        colour = "grey60", alpha = 0.3) +
    ggplot2::geom_point(data = new,
                        ggplot2::aes(colour = .data$within_domain),
                        size = 2.5) +
    ggplot2::geom_text(data = labels,
                       ggplot2::aes(label = .data$label, hjust = .data$hjust)) +
    ggplot2::scale_size_continuous(range = c(0.3, 3),
                                   name = "Training NPP") +
    scale_domain() +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(title = "brGDGT methylation groups") +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom",
                   plot.margin = ggplot2::margin(10, 30, 10, 30))
}

plot_similarity_map <- function(object, n_sites = 10) {
  if (!requireNamespace("maps", quietly = TRUE)) {
    stop("Package 'maps' is needed for the similarity map. ",
         "Install it with install.packages(\"maps\")", call. = FALSE)
  }
  # training sites ranked by similarity: 1 / (1 + mean distance to the
  # samples for which they are nearest neighbours)
  nb <- object$neighbours
  sim <- tapply(nb$distance, nb$ref_index, function(d) 1 / (1 + mean(d)))
  top <- as.integer(names(sort(sim, decreasing = TRUE)))
  top <- top[seq_len(min(n_sites, length(top)))]

  coords <- object$reference$coordinates
  all_sites <- data.frame(longitude = coords$Longitude,
                          latitude = coords$Latitude)
  top_sites <- data.frame(longitude = coords$Longitude[top],
                          latitude = coords$Latitude[top],
                          similarity = as.numeric(sim[as.character(top)]))
  top_sites <- top_sites[stats::complete.cases(top_sites), ]
  world <- ggplot2::map_data("world")

  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$longitude,
                                         y = .data$latitude)) +
    ggplot2::geom_polygon(data = world,
                          ggplot2::aes(x = .data$long, y = .data$lat,
                                       group = .data$group),
                          fill = "grey90", colour = "white", linewidth = 0.2) +
    ggplot2::geom_point(data = all_sites, colour = "black", size = 0.5,
                        alpha = 0.3) +
    ggplot2::geom_point(data = top_sites,
                        ggplot2::aes(size = .data$similarity),
                        colour = "darkorange", alpha = 0.7) +
    ggplot2::scale_size_continuous(range = c(3, 8), guide = "none") +
    ggplot2::coord_fixed(1.3) +
    ggplot2::labs(title = "Training sites with the most similar brGDGTs",
                  subtitle = "Euclidean distance in CLR space",
                  x = "Longitude", y = "Latitude") +
    ggplot2::theme_minimal()
}

#' Predict potential NPP from climate
#'
#' @description Predicts potential Net Primary Productivity (NPP) from mean
#'   annual temperature and precipitation.
#'
#'   `npp_miami()` uses the Miami model (Lieth 1975): NPP is the lower of a
#'   temperature-limited and a precipitation-limited estimate.
#'
#'   `npp_thornthwaite()` uses the Thornthwaite Memorial model (Lieth 1972),
#'   which predicts NPP from actual evapotranspiration (AET), here estimated
#'   from temperature and precipitation with Turc's formula. NPP is set to 0
#'   where AET is below 20 mm.
#'
#' @param mat Mean annual temperature (degrees C).
#' @param map Mean annual precipitation (mm), same length as `mat`.
#' @param id Optional sample identifiers; defaults to `seq_along(mat)`.
#'
#' @return An `npp_climate` data frame with columns `id`, `mat`, `map`, `npp`
#'   (gC m-2 yr-1) and model-specific columns: `npp_t`, `npp_p` and
#'   `limiting` for the Miami model, `aet` for the Thornthwaite model.
#'
#' @references Lieth, H. (1975) Modeling the primary productivity of the
#'   world. In: Primary Productivity of the Biosphere, 237-263.
#'   doi:10.1007/978-3-642-80913-2_12
#'
#'   Lieth, H. & Box, E. (1972) Evapotranspiration and primary productivity:
#'   C.W. Thornthwaite memorial model. Publications in Climatology 25, 37-46.
#' @export
#'
#' @examples
#' npp_miami(mat = c(5, 10, 20), map = c(400, 800, 1500))
#' npp_thornthwaite(mat = c(5, 10, 20), map = c(400, 800, 1500))
npp_miami <- function(mat, map, id = NULL) {
  id <- check_climate(mat, map, id)
  # Two empirical curves fitted by Lieth to global NPP data, both saturating
  # at 3000 gC m-2 yr-1:
  # a logistic (S-shaped) curve in temperature...
  npp_t <- 3000 / (1 + exp(1.315 - 0.119 * mat))
  # ...and a saturating exponential in precipitation
  npp_p <- 3000 * (1 - exp(-0.000664 * map))
  # Liebig's law of the minimum: productivity is set by whichever factor is
  # more limiting. Both curves are returned, together with which one was
  # limiting, as that is often of interest in itself.
  out <- data.frame(
    id = id, mat = mat, map = map,
    npp = pmin(npp_t, npp_p),
    npp_t = npp_t,
    npp_p = npp_p,
    limiting = ifelse(npp_t <= npp_p, "temperature", "precipitation")
  )
  new_npp_climate(out, "Miami")
}

#' @rdname npp_miami
#' @export
npp_thornthwaite <- function(mat, map, id = NULL) {
  id <- check_climate(mat, map, id)
  # Turc's formula for actual evapotranspiration (AET, mm/yr):
  #   L = 300 + 25 T + 0.05 T^3   (evaporating power of the air, mm/yr)
  #   AET = P / sqrt(0.9 + (P / L)^2)
  # L becomes negative below about -10 C, where there is effectively no
  # evaporation, so it is floored at 0.
  evap_max <- pmax(300 + 25 * mat + 0.05 * mat^3, 0)
  # When P is small relative to L (P < 0.316 L), Turc's formula would give
  # AET > P, i.e. more water evaporating than falls. AET is capped at P there.
  # The two branches meet at P = 0.316 L, since sqrt(0.9 + 0.316^2) = 1, so
  # the curve is continuous. Using <= sends P = L = 0 (no rain in a very cold
  # climate) to the first branch, giving AET = 0 rather than 0 / 0 = NaN.
  aet <- ifelse(map <= 0.316 * evap_max,
                map,
                map / sqrt(0.9 + (map / evap_max)^2))
  # Lieth's Thornthwaite Memorial curve: NPP saturates at 3000 gC m-2 yr-1 as
  # AET increases. The curve crosses zero at AET = 20 mm and would go negative
  # below that, which is meaningless, so it is floored at 0.
  out <- data.frame(
    id = id, mat = mat, map = map,
    npp = pmax(3000 * (1 - exp(-0.0009695 * (aet - 20))), 0),
    aet = aet
  )
  new_npp_climate(out, "Thornthwaite")
}

# Validate the climate inputs shared by both models, and return the ids to
# use (sequential numbers if the user gave none)
check_climate <- function(mat, map, id) {
  if (!is.numeric(mat) || !is.numeric(map)) {
    stop("`mat` and `map` must be numeric", call. = FALSE)
  }
  if (length(mat) != length(map)) {
    stop("`mat` and `map` must have the same length", call. = FALSE)
  }
  if (is.null(id)) id <- seq_along(mat)
  if (length(id) != length(mat)) {
    stop("`id` must have the same length as `mat`", call. = FALSE)
  }
  id
}

# Constructor for the climate NPP class. The model name is kept as an
# attribute so that print(), autoplot() and hanpp() can report it.
new_npp_climate <- function(df, model) {
  structure(df, class = c("npp_climate", "data.frame"), model = model)
}

#' @export
print.npp_climate <- function(x, ...) {
  cat("Potential NPP (gC m-2 yr-1) from the", attr(x, "model"),
      "model for", nrow(x), "samples\n\n")
  print(as.data.frame(x), ...)
  invisible(x)
}

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
  if (is.null(xlab)) xlab <- if (is.null(x)) "Sample" else deparse(substitute(x))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$npp)) +
    ggplot2::geom_line(colour = "grey50")
  # only the Miami model reports a limiting factor, so colour by it if present
  if ("limiting" %in% names(df)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$limiting)) +
      ggplot2::labs(colour = "Limiting factor")
  } else {
    p <- p + ggplot2::geom_point()
  }
  p + ggplot2::labs(x = xlab, y = "NPP (gC m-2 yr-1)",
                    title = paste("Potential NPP,", attr(object, "model"),
                                  "model")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

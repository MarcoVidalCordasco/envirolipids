#' Predict potential NPP from climate
#'
#' @description Predicts potential Net Primary Productivity (NPP) from mean
#'   annual temperature and total annual precipitation.
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
#' @param map Total annual precipitation (mm per year), as a long-term average
#'   of yearly totals; not mean monthly precipitation. Same length as `mat`.
#' @param id Optional sample identifiers; defaults to `seq_along(mat)`.
#' @param units Units of the returned NPP. `"carbon"` (the default) gives
#'   grams of carbon (gC m-2 yr-1), the same units as [npp_brgdgt()], so that
#'   the two can be compared. `"dry_matter"` gives grams of dry matter
#'   (g m-2 yr-1), the units of Lieth's original models. Carbon is taken to
#'   be 45% of dry matter.
#'
#' @return An `npp_climate` data frame with columns `id`, `mat`, `map`, `npp`
#'   and model-specific columns: `npp_t`, `npp_p` and `limiting` for the Miami
#'   model, `aet` for the Thornthwaite model. All NPP columns are in the
#'   requested `units`.
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
npp_miami <- function(mat, map, id = NULL,
                      units = c("carbon", "dry_matter")) {
  id <- check_climate(mat, map, id)
  units <- match.arg(units)
  # Lieth's models predict NPP as grams of dry matter (g m-2 yr-1); by default
  # this is converted to grams of carbon (see dm_to_units()) to match
  # brGDGT-based NPP.
  #
  # Two empirical curves fitted by Lieth to global NPP data, both saturating
  # at 3000 g m-2 yr-1:
  # a logistic (S-shaped) curve in temperature...
  npp_t <- 3000 / (1 + exp(1.315 - 0.119 * mat))
  # ...and a saturating exponential in precipitation
  npp_p <- 3000 * (1 - exp(-0.000664 * map))
  # Liebig's law of the minimum: productivity is set by whichever factor is
  # more limiting. Both curves are returned, together with which one was
  # limiting, as that is often of interest in itself.
  out <- data.frame(
    id = id, mat = mat, map = map,
    npp = dm_to_units(pmin(npp_t, npp_p), units),
    npp_t = dm_to_units(npp_t, units),
    npp_p = dm_to_units(npp_p, units),
    limiting = ifelse(npp_t <= npp_p, "temperature", "precipitation")
  )
  new_npp_climate(out, "Miami", units)
}

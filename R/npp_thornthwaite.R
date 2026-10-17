#' @rdname npp_miami
#' @export
npp_thornthwaite <- function(mat, map, id = NULL,
                             units = c("carbon", "dry_matter")) {
  id <- check_climate(mat, map, id)
  units <- match.arg(units)
  # Turc's formula for actual evapotranspiration (AET, mm/yr):
  # nolint start: commented_code_linter.
  #   L = 300 + 25 T + 0.05 T^3   (evaporating power of the air, mm/yr)
  #   AET = P / sqrt(0.9 + (P / L)^2)
  # nolint end
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
  # Lieth's Thornthwaite Memorial curve: NPP saturates at 3000 g m-2 yr-1 as
  # AET increases. The curve crosses zero at AET = 20 mm and would go negative
  # below that, which is meaningless, so it is floored at 0.
  out <- data.frame(
    id = id, mat = mat, map = map,
    npp = dm_to_units(pmax(3000 * (1 - exp(-0.0009695 * (aet - 20))), 0),
                      units),
    aet = aet
  )
  new_npp_climate(out, "Thornthwaite", units)
}

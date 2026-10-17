# Internal helpers shared by npp_miami() and npp_thornthwaite(), and the
# print method for the npp_climate class they return.

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

# Convert NPP from grams of dry matter, as given by Lieth's models, to the
# requested units. Carbon is taken to be 45% of dry matter, a standard
# conversion factor for plant biomass.
dm_to_units <- function(npp, units) {
  if (units == "carbon") npp * 0.45 else npp
}

# Constructor for the climate NPP class. The model name and the units label
# are kept as attributes so that print() and autoplot() can report them.
new_npp_climate <- function(df, model, units) {
  structure(df, class = c("npp_climate", "data.frame"), model = model,
            units = if (units == "carbon") "gC m-2 yr-1" else "g m-2 yr-1")
}

#' @export
print.npp_climate <- function(x, ...) {
  cat("Potential NPP (", attr(x, "units"), ") from the ", attr(x, "model"),
      " model for ", nrow(x), " samples\n\n", sep = "")
  print(as.data.frame(x), ...)
  invisible(x)
}

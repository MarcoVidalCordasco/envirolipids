#' Predict air-temperature metrics from brGDGTs
#'
#' @description
#' Select the sample material first, then choose a published temperature
#' calibration appropriate to that material. Published equations estimate
#' air-temperature metrics.
#'
#' Available sample types are `"soils"`, `"mineral_soils"`, and
#' `"lake_sediments"`. The first category offers the general soil
#' calibrations (Peterse and De Jonge). `"mineral_soils"` additionally offers
#' the Naafs calibrations developed for mineral soils. `"lake_sediments"`
#' offers lake-sediment calibrations (Russell, Raberg and Zhu).
#'
#' MAT = mean annual air temperature (deg C); MAF = mean temperature of months
#' above freezing (deg C); GDD0 = growing degree days above 0 C (deg C days/year).
#'
#' Call `temp_brgdgt_methods(sample_type)` to see the methods available for
#' your material, then pass one of those names as `method`.
#'
#' @param data A data frame with the 15 brGDGT compounds as fractions or
#'   percentages.
#' @param sample_type One sample-material category: `"soils"`,
#'   `"mineral_soils"`, or `"lake_sediments"`. Spaces and singular forms such
#'   as `"mineral soils"` and `"lake sediment"` are also accepted. Each call
#'   assumes the samples share the selected material category.
#' @param method One of the method names returned by
#'   `temp_brgdgt_methods(sample_type)`. Required; no calibration is selected
#'   silently. The special `"random_forest"` method remains available if a
#'   fitted temperature model is supplied in `model`.
#' @param ci For `method = "random_forest"`, calculate tree-quantile bounds.
#'   Linear calibration RMSEs are not confidence intervals.
#' @param id_col Name of the sample identifier column, default `"ID"`.
#' @param model Fitted random forest model, only for `method = "random_forest"`.
#'
#' @return A `brgdgt_prediction` data frame giving the estimate, target
#'   (MAT / MAF / GDD0), units, publication, equation, calibration material,
#'   selected sample material, and published RMSE. The restricted Naafs
#'   equations set estimates to `NA` for samples with IR6me >= 0.5.
#'   Other environmental limitations (for example freshwater-only lake
#'   calibrations) remain the user's responsibility.
#'
#' @references
#' Peterse et al. (2012); De Jonge et al. (2014);
#' Naafs et al. (2017), Organic Geochemistry 106, 48-56;
#' Russell et al. (2018); Raberg et al. (2021), Biogeosciences 18, 3579-3603;
#' Zhu et al. (2025), Quaternary Science Reviews 357, 109319.
#'
#' @examples
#' data(example_data)
#' temp_brgdgt_methods("mineral_soils")
#' temp_brgdgt_methods("lake_sediments")
#' soil_mat <- temp_brgdgt(example_data, sample_type = "mineral_soils",
#'                         method = "naafs2017_global")
#' lake_maf <- temp_brgdgt(example_data, sample_type = "lake_sediments",
#'                         method = "zhu2025_1")
#' @export
temp_brgdgt <- function(data, sample_type, method = NULL, ci = FALSE,
                        id_col = "ID", model = NULL) {
  if (missing(sample_type)) {
    stop("Choose `sample_type` first: soils, mineral_soils, or ",
         "lake_sediments. See temp_brgdgt_methods().", call. = FALSE)
  }
  sample_type <- .temp_normalize_sample_type(sample_type)
  available <- .temp_methods_for_sample(sample_type)
  
  if (is.null(method)) {
    stop("Please choose `method` for sample_type = '", sample_type,
         "'. Available methods: ", paste(available$method, collapse = ", "),
         ". You can also call temp_brgdgt_methods('", sample_type, "').",
         call. = FALSE)
  }
  if (!is.character(method) || length(method) != 1L ||
      is.na(method) || !nzchar(method)) {
    stop("`method` must be one non-missing method name", call. = FALSE)
  }
  if (!is.logical(ci) || length(ci) != 1L || is.na(ci)) {
    stop("`ci` must be TRUE or FALSE", call. = FALSE)
  }
  
  # Backward-compatible optional RF route. This is not one of the published
  # equations in the method catalogue; the original package did not ship a
  # temperature model. Translate the new material categories to the original
  # predictor categories expected by predict_rf_brgdgt().
  if (identical(method, "random_forest")) {
    if (is.null(model)) {
      stop("A fitted temperature random-forest `model` must be supplied; ",
           "the package does not include one.", call. = FALSE)
    }
    rf_type <- if (sample_type == "lake_sediments") {
      "Lacustrine Sediment"
    } else {
      "Soil"
    }
    out <- predict_rf_brgdgt(data, rf_type, model = model, ci = ci,
                             id_col = id_col, response = "MAT",
                             units = "\u00b0C")
    out$material_group <- sample_type
    return(out)
  }
  
  if (!(method %in% available$method)) {
    stop("Method '", method, "' is not available for sample_type = '",
         sample_type, "'. Available methods: ",
         paste(available$method, collapse = ", "),
         ". See temp_brgdgt_methods('", sample_type, "').",
         call. = FALSE)
  }
  if (ci) {
    warning("`ci` applies only to random forests; published calibration ",
            "RMSE is reported instead (not a confidence interval).",
            call. = FALSE)
  }
  
  fractions <- brgdgt_fractions(data)
  out <- predict_temp_calibration(fractions, method)
  out$id <- get_ids(data, id_col)
  out$sample_type <- rep(sample_type, nrow(out))
  # All shown published methods are compatible with the chosen broad
  # material category. This does not imply all regional or chemical
  # applicability assumptions of a specific calibration are met.
  out$environment_match <- rep(TRUE, nrow(out))
  out <- out[, c("id", "sample_type", setdiff(names(out),
                                              c("id", "sample_type"))), drop = FALSE]
  rownames(out) <- NULL
  
  if (any(!out$applicable)) {
    warning(sum(!out$applicable),
            " sample(s) fail the IR6me < 0.5 restriction; ",
            "their estimates were set to NA.", call. = FALSE)
  }
  
  meta <- .temp_catalogue()
  ref <- meta[meta$method == method, , drop = FALSE]
  structure(out, class = c("brgdgt_prediction", "data.frame"),
            response = ref$target, units = ref$units,
            method = method, sample_type = sample_type)
}

# Internal: convert friendly sample-type labels to canonical codes.
.temp_normalize_sample_type <- function(sample_type) {
  allowed <- c("soils", "mineral_soils", "lake_sediments")
  if (!is.character(sample_type) || length(sample_type) != 1L ||
      is.na(sample_type) || !nzchar(trimws(sample_type))) {
    stop("`sample_type` must be one of: ", paste(allowed, collapse = ", "),
         call. = FALSE)
  }
  label <- tolower(trimws(sample_type))
  label <- gsub("[[:space:]-]+", "_", label)
  aliases <- c(soil = "soils", soils = "soils",
               mineral_soil = "mineral_soils",
               mineral_soils = "mineral_soils",
               lake_sediment = "lake_sediments",
               lake_sediments = "lake_sediments",
               lacustrine_sediment = "lake_sediments",
               lacustrine_sediments = "lake_sediments")
  if (label %in% names(aliases)) label <- unname(aliases[[label]])
  if (!(label %in% allowed)) {
    stop("Unknown `sample_type`: ", sample_type,
         ". Choose one of: ", paste(allowed, collapse = ", "),
         call. = FALSE)
  }
  label
}

# Internal: the general soil calibrations also apply to mineral soils.
.temp_methods_for_sample <- function(sample_type) {
  sample_type <- .temp_normalize_sample_type(sample_type)
  allowed_groups <- if (identical(sample_type, "mineral_soils")) {
    c("soils", "mineral_soils")
  } else {
    sample_type
  }
  catalogue <- .temp_catalogue()
  out <- catalogue[catalogue$sample_group %in% allowed_groups, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' List brGDGT temperature calibrations for a selected sample material
#'
#' @description
#' Discover temperature calibrations by sample material, prior to calling
#' `temp_brgdgt()`. The `"soils"` category lists general soil equations;
#' `"mineral_soils"` includes these plus Naafs et al. (2017); and
#' `"lake_sediments"` lists lake calibrations. Without `sample_type`,
#' all 18 implemented literature equations are returned.
#'
#' @param sample_type Optional material: `"soils"`, `"mineral_soils"`, or
#'   `"lake_sediments"`. `NULL` (the default) returns all published methods.
#' @return A data frame of available method names and scientific metadata:
#'   target temperature (MAT, MAF, GDD0), units, calibration material,
#'   publication, equation number, reported RMSE, and limitations.
#' @examples
#' temp_brgdgt_methods("soils")
#' temp_brgdgt_methods("mineral_soils")
#' temp_brgdgt_methods("lake_sediments")
#' @export
temp_brgdgt_methods <- function(sample_type = NULL) {
  if (is.null(sample_type)) return(.temp_catalogue())
  .temp_methods_for_sample(sample_type)
}

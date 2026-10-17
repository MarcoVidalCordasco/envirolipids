# Build data/modern_example_data.rda: a small set of modern samples taken from
# the training dataset (inst/extdata/Dataset_.xlsx), with present-day climate
# from WorldClim 2.1 (1970-2000 averages, 10 arc-minute resolution).
# Run from the package root. Requires readxl and terra, which are not
# package dependencies, and downloads about 50 MB from WorldClim the first
# time it is run.
#
# Usage: Rscript data-raw/make_modern_example_data.R [worldclim_dir]
# where worldclim_dir holds (or will receive) wc2.1_10m_bio.zip.

brgdgt_cols <- c("fIa", "fIb", "fIc", "fIIa", "fIIa_", "fIIb", "fIIb_",
                 "fIIc", "fIIc_", "fIIIa", "fIIIa_", "fIIIb", "fIIIb_",
                 "fIIIc", "fIIIc_")

args <- commandArgs(trailingOnly = TRUE)
wc_dir <- if (length(args) > 0) args[1] else file.path(tempdir(), "worldclim")
dir.create(wc_dir, showWarnings = FALSE, recursive = TRUE)

# --- WorldClim layers --------------------------------------------------------
# bio_1 = mean annual temperature (degrees C); bio_12 = total annual
# precipitation (mm). In WorldClim 2.1 these are stored in real units (unlike
# version 1.4, where temperatures were multiplied by 10).
wc_zip <- file.path(wc_dir, "wc2.1_10m_bio.zip")
if (!file.exists(wc_zip)) {
  utils::download.file(
    "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip",
    wc_zip, mode = "wb"
  )
}
layers <- c(mat = "wc2.1_10m_bio_1.tif", map = "wc2.1_10m_bio_12.tif")
utils::unzip(wc_zip, files = layers, exdir = wc_dir)
climate <- terra::rast(file.path(wc_dir, layers))
names(climate) <- names(layers)

# --- Candidate samples -------------------------------------------------------
raw_data <- as.data.frame(
  readxl::read_excel("inst/extdata/Dataset_.xlsx", sheet = "Dataset")
)
raw_data$Site <- trimws(raw_data$Site)

# Only samples that could have been used to train the NPP model (those with
# observed NPP) and with valid coordinates. This also drops the 104 samples
# whose latitude and longitude are swapped in the spreadsheet (latitude > 90).
ok <- !is.na(raw_data$NPP) &
  !is.na(raw_data$Latitude) & abs(raw_data$Latitude) <= 90 &
  !is.na(raw_data$Longitude) & abs(raw_data$Longitude) <= 180
candidates <- raw_data[ok, ]

# Climate at each sample location. Points on water or on coasts can fall on
# cells without data, so keep only samples with both values.
xy <- as.matrix(candidates[, c("Longitude", "Latitude")])
clim <- terra::extract(climate, xy)
candidates$MAT <- round(clim$mat, 1)
candidates$MAP <- round(clim$map)
candidates <- candidates[!is.na(candidates$MAT) & !is.na(candidates$MAP), ]

# One sample per location, so that the example covers distinct sites
duplicated_site <- duplicated(candidates[, c("Latitude", "Longitude")])
candidates <- candidates[!duplicated_site, ]

# --- Selection ---------------------------------------------------------------
# A fixed mix of the three most common sample types. Within each type, samples
# from different studies (Reference) come first, which spreads the selection
# across studies and so across regions. Types with few studies (e.g. peat) are
# then topped up with further samples.
n_per_type <- c("Soil" = 10, "Lacustrine Sediment" = 6, "Peat" = 4)
set.seed(2026)
selected <- do.call(rbind, lapply(names(n_per_type), function(type) {
  pool <- candidates[candidates$Sampletype == type, ]
  pool <- pool[sample(nrow(pool)), ]
  # order() is stable, so this keeps the random order within each group
  pool <- pool[order(duplicated(pool$Reference)), ]
  pool[seq_len(min(n_per_type[[type]], nrow(pool))), ]
}))

# Accented characters in the spreadsheet were corrupted to "ï¿½"
# at some point, and the original letters cannot be recovered automatically.
# Fix the selected site names by hand (ASCII only, as CRAN requires for data),
# and stop if a new selection brings in a garbled name not covered here.
site_fixes <- c("Fjallsjï¿½kull" = "Fjallsjokull")
to_fix <- selected$Site %in% names(site_fixes)
selected$Site[to_fix] <- site_fixes[selected$Site[to_fix]]
stopifnot(!grepl("[^ -~]", selected$Site))

modern_example_data <- data.frame(
  ID = selected$Site,
  Sampletype = selected$Sampletype,
  Latitude = selected$Latitude,
  Longitude = selected$Longitude,
  # keep percentages as in the spreadsheet, to show that the package accepts
  # them as well as fractions
  selected[, brgdgt_cols],
  NPP = selected$NPP,
  MAT = selected$MAT,
  MAP = selected$MAP,
  row.names = NULL
)

save(modern_example_data, file = "data/modern_example_data.rda",
     compress = "bzip2")

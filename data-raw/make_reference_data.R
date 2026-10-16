# Build inst/extdata/reference_data.rds from the training dataset
# (inst/extdata/Dataset_.xlsx). Run from the package root.
# Requires readxl, which is not a package dependency.

brgdgt_cols <- c("fIa", "fIb", "fIc", "fIIa", "fIIa_", "fIIb", "fIIb_",
                 "fIIc", "fIIc_", "fIIIa", "fIIIa_", "fIIIb", "fIIIb_",
                 "fIIIc", "fIIIc_")

raw_data <- as.data.frame(
  readxl::read_excel("inst/extdata/Dataset_.xlsx", sheet = "Dataset")
)

# Training abundances are percentages: convert to fractions, then apply the
# CLR transform with the same 1e-6 zero offset as brgdgt_clr() in
# R/brgdgt_utils.R. The two must stay identical, or new samples and training
# data would be transformed differently.
fractions <- raw_data[, brgdgt_cols] / 100
log_comp <- log(fractions + 1e-6)
clr_values <- as.data.frame(log_comp - rowMeans(log_comp))

# Store both fractions (for methylation summaries and plots) and CLR values
# (for distances), so they are not recomputed on every call. Site names are
# trimmed, as some have trailing spaces in the spreadsheet. NPP is missing (NA)
# for some samples.
reference_data <- list(
  fractions = fractions,
  clr_values = clr_values,
  npp = raw_data$NPP,
  sample_types = raw_data$Sampletype,
  sites = trimws(raw_data$Site),
  coordinates = raw_data[, c("Longitude", "Latitude")]
)

saveRDS(reference_data, "inst/extdata/reference_data.rds")

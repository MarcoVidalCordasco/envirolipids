#' @keywords internal
"_PACKAGE"

#' @importFrom ggplot2 autoplot .data
#' @importFrom randomForest randomForest
NULL

# Constants shared across the package ----------------------------------------

# The 15 brGDGT compounds used by the models. All functions select the
# compounds from user data by these names (never by column position), so the
# order of columns in the user's data does not matter
BRGDGT_COLS <- c("fIa", "fIb", "fIc", "fIIa", "fIIa_", "fIIb", "fIIb_",
                 "fIIc", "fIIc_", "fIIIa", "fIIIa_", "fIIIb", "fIIIb_",
                 "fIIIc", "fIIIc_")

# Compounds grouped by the number of methyl groups on the alkyl chains:
# I = tetramethylated, II = pentamethylated, III = hexamethylated
# (a trailing underscore denotes the 6-methyl isomer, e.g. IIa')
METHYLATION_GROUPS <- list(
  tetra = c("fIa", "fIb", "fIc"),
  penta = c("fIIa", "fIIa_", "fIIb", "fIIb_", "fIIc", "fIIc_"),
  hexa = c("fIIIa", "fIIIa_", "fIIIb", "fIIIb_", "fIIIc", "fIIIc_")
)

# Sample types present in the training data, i.e. the values users can pass
# as `sample_type`
SAMPLE_TYPES <- c("Lacustrine Sediment", "Lacustrine SPM",
                  "Low DO Lacustrine SPM", "Peat",
                  "Riverine Sediment and SPM", "Soil",
                  "Lacustrine Meso/Microcosm")

# Sample types pooled into a single class when training the models
AQUATIC_SPM_TYPES <- c("Lacustrine SPM", "Low DO Lacustrine SPM",
                       "Riverine Sediment and SPM", "Lacustrine Meso/Microcosm")

# Colours for samples within (green) and outside (red) the training domain,
# shared by all plots so the legend means the same thing everywhere
DOMAIN_COLOURS <- c("TRUE" = "#2E8B57", "FALSE" = "#DC143C")

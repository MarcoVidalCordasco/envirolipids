#' Example paleo brGDGT dataset from a sediment core
#'
#' A dataset containing brGDGT measurements from 50 natural archive samples
#' to demonstrate the package functionality on paleo data. See
#' [modern_example_data] for modern samples. The data includes fractional
#' abundances of 15 brGDGT compounds and sample ages (years BP).
#'
#' @format A data frame with 50 rows and 17 variables (sample ID, 15 brGDGTs
#'   and Age):
#' \describe{
#'   \item{ID}{Sample identifier (numeric)}
#'   \item{fIIIa}{brGDGT IIIa fractional abundance}
#'   \item{fIIIa_}{brGDGT IIIa' fractional abundance}
#'   \item{fIIIb}{brGDGT IIIb fractional abundance}
#'   \item{fIIIb_}{brGDGT IIIb' fractional abundance}
#'   \item{fIIIc}{brGDGT IIIc fractional abundance}
#'   \item{fIIIc_}{brGDGT IIIc' fractional abundance}
#'   \item{fIIa}{brGDGT IIa fractional abundance}
#'   \item{fIIa_}{brGDGT IIa' fractional abundance}
#'   \item{fIIb}{brGDGT IIb fractional abundance}
#'   \item{fIIb_}{brGDGT IIb' fractional abundance}
#'   \item{fIIc}{brGDGT IIc fractional abundance}
#'   \item{fIIc_}{brGDGT IIc' fractional abundance}
#'   \item{fIa}{brGDGT Ia fractional abundance}
#'   \item{fIb}{brGDGT Ib fractional abundance}
#'   \item{fIc}{brGDGT Ic fractional abundance}
#'   \item{Age}{Sample age (years BP)}
#' }
#'
#' @details
#' The dataset contains 50 samples from a sediment core. All brGDGT abundances
#' are expressed as fractions (summing to 1 per sample).
#'
#' @source Natural archive samples. For demonstration purposes only.
#'
#' @examples
#' \dontrun{
#' # Load the example data
#' data(example_data)
#'
#' # View structure
#' str(example_data)
#'
#' # Check first few rows
#' head(example_data)
#'
#' # Calculate sum of brGDGTs (should be 1)
#' brGDGT_cols <- c("fIa", "fIb", "fIc", "fIIa", "fIIa_", "fIIb", "fIIb_",
#'                  "fIIc", "fIIc_", "fIIIa", "fIIIa_", "fIIIb", "fIIIb_",
#'                  "fIIIc", "fIIIc_")
#' rowSums(example_data[, brGDGT_cols])
#' }
#'
"example_data"

#' Example modern brGDGT dataset
#'
#' Twenty modern samples (10 soils, 6 lake sediments and 4 peats) from sites
#' around the world, to demonstrate the package on modern data. Each sample
#' comes with its observed NPP and present-day climate.
#'
#' @format A data frame with 20 rows and 22 variables:
#' \describe{
#'   \item{ID}{Site name}
#'   \item{Sampletype}{Sample type, as accepted by [npp_brgdgt()]}
#'   \item{Latitude, Longitude}{Site coordinates (decimal degrees)}
#'   \item{fIa, fIb, fIc, fIIa, fIIa_, fIIb, fIIb_, fIIc, fIIc_, fIIIa, fIIIa_,
#'     fIIIb, fIIIb_, fIIIc, fIIIc_}{brGDGT abundances, as percentages (summing
#'     to 100 per sample)}
#'   \item{NPP}{Observed NPP at the site (gC m-2 yr-1)}
#'   \item{MAT}{Mean annual temperature (degrees C), from WorldClim 2.1}
#'   \item{MAP}{Total annual precipitation (mm per year), from WorldClim 2.1}
#' }
#'
#' @details The samples are a subset of the data used to train the NPP model,
#'   selected to cover different studies and regions. Because the model has
#'   seen these samples, its predictions for them are better than can be
#'   expected for genuinely new samples. Climate is the WorldClim 2.1 average
#'   for 1970-2000 at 10 arc-minute resolution, extracted at each site's
#'   coordinates. The dataset is built by `data-raw/make_modern_example_data.R`.
#'
#' @source Training compilation of the NPP model; climate from Fick, S.E. &
#'   Hijmans, R.J. (2017) WorldClim 2: new 1-km spatial resolution climate
#'   surfaces for global land areas. International Journal of Climatology 37,
#'   4302-4315. doi:10.1002/joc.5086
#'
#' @examples
#' data(modern_example_data)
#' head(modern_example_data[, c("ID", "Sampletype", "NPP", "MAT", "MAP")])
"modern_example_data"

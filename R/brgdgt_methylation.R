#' Summarise brGDGTs by methylation group
#'
#' @description Sums the fractional abundances of tetra-, penta- and
#'   hexamethylated brGDGTs for each sample.
#'
#' @param data A data frame with the 15 brGDGT columns (fIa, fIb, fIc, fIIa,
#'   fIIa_, fIIb, fIIb_, fIIc, fIIc_, fIIIa, fIIIa_, fIIIb, fIIIb_, fIIIc,
#'   fIIIc_), as fractions or percentages.
#' @param id_col Name of the column holding sample identifiers. If absent,
#'   samples are numbered by row.
#'
#' @return A data frame with columns `id`, `tetra`, `penta` and `hexa`; each
#'   row sums to 1.
#' @export
#'
#' @examples
#' data(example_data)
#' head(brgdgt_methylation(example_data))
brgdgt_methylation <- function(data, id_col = "ID") {
  fractions <- brgdgt_fractions(data)
  data.frame(id = get_ids(data, id_col), methylation_groups(fractions))
}

################################################################################
#' Count taxa and sequences by congruence of two taxonomic assignments
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Classifies each taxon of a phyloseq object into a congruence category by
#' comparing two `tax_table` columns (typically the same taxonomic rank assigned
#' by two different databases, algorithms or reference sets), and counts both the
#' number of taxa (ASV/OTU) and the number of sequences (reads) falling in each
#' category. The five categories are: `both_equal` (identical non-missing
#' assignment), `both_na` (missing in both), `only_<rank_1>` (assigned only by
#' the first column), `only_<rank_2>` (assigned only by the second column), and
#' `different` (both assigned but disagreeing). Empty strings and the literal
#' values `"NA"` / `"NA_NA"` are treated as missing.
#'
#' @param physeq (required) A [phyloseq::phyloseq-class()] object.
#' @param rank_1 (character, required) Name of the `tax_table` column for the
#'   first assignment (e.g. `"Genus"` or `"Genus_SILVA"`).
#' @param rank_2 (character, required) Name of the `tax_table` column for the
#'   second assignment (e.g. `"Genus__eukaryome_Glomero"`).
#'
#' @returns A [tibble::tibble()] with one row per category and the columns
#'   `category`, `n_asv` (number of taxa), `n_seq` (summed sequence counts),
#'   `pct_asv` and `pct_seq` (percentages, rounded to two decimals). The
#'   `only_<rank_1>` and `only_<rank_2>` category names embed the actual column
#'   names passed as arguments.
#' @export
#' @author Adrien Taudière
#'
#' @seealso [plot_congruence_counts()], [compare_taxo_db()]
#'
#' @examples
#' \donttest{
#' count_taxo_congruence(Glom_otu, "Genus", "Genus__eukaryome_Glomero")
#' }
count_taxo_congruence <- function(physeq, rank_1, rank_2) {
  verify_pq(physeq)

  clean_na <- function(x) {
    x[x == "" | x == "NA_NA" | x == "NA"] <- NA
    x
  }

  v1 <- clean_na(as.character(physeq@tax_table[, rank_1]))
  v2 <- clean_na(as.character(physeq@tax_table[, rank_2]))

  na1 <- is.na(v1)
  na2 <- is.na(v2)

  only_1 <- paste0("only_", rank_1)
  only_2 <- paste0("only_", rank_2)

  category <- ifelse(
    na1 & na2,
    "both_na",
    ifelse(
      !na1 & na2,
      only_1,
      ifelse(
        na1 & !na2,
        only_2,
        ifelse(v1 == v2, "both_equal", "different")
      )
    )
  )

  seq_sums <- taxa_sums(physeq)

  cats <- c("both_equal", "both_na", only_1, only_2, "different")
  result <- data.frame(
    category = cats,
    n_asv = vapply(cats, function(c) sum(category == c), numeric(1)),
    n_seq = vapply(cats, function(c) sum(seq_sums[category == c]), numeric(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  ) |>
    tibble::tibble() |>
    dplyr::mutate(
      pct_asv = round(100 * n_asv / sum(n_asv), 2),
      pct_seq = round(100 * n_seq / sum(n_seq), 2)
    )
  result
}
################################################################################

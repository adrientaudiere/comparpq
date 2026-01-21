################################################################################
#' Contingency table of two taxonomic ranks
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Creates a cross-tabulation (contingency table) comparing two taxonomic ranks
#' from a phyloseq object. Useful for comparing taxonomic assignments from
#' different databases, algorithms, or taxonomic levels.
#'
#' @inheritParams tc_points_matrix
#' @param rank_1 (character, default "Family") The name of the first taxonomic
#'   rank (column in tax_table slot) for the cross-tabulation rows.
#' @param rank_2 (character, default "Class") The name of the second taxonomic
#'   rank (column in tax_table slot) for the cross-tabulation columns.
#' @param ... Additional arguments passed to [gtsummary::tbl_cross()].
#'
#' @returns A gtsummary tbl_cross object displaying the cross-tabulation of
#'   the two taxonomic ranks.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' tc_df_pq(data_fungi_mini)
#' tc_df_pq(data_fungi_mini, rank_1 = "Order", rank_2 = "Family")
#'
#' \dontrun{
#' # Compare taxonomic assignments from different methods
#' ref_fasta <- system.file("extdata", "mini_UNITE_fungi.fasta.gz",
#'   package = "MiscMetabar", mustWork = TRUE
#' )
#' data_fungi_mini2 <- data_fungi_mini |>
#'   add_new_taxonomy_pq(ref_fasta, suffix = "_sintax", method = "sintax") |>
#'   add_new_taxonomy_pq(ref_fasta, suffix = "_lca", method = "lca")
#'
#' tc_df_pq(data_fungi_mini2, rank_1 = "Class_lca", rank_2 = "Class_sintax")
#' }
tc_df_pq <- function(physeq, rank_1 = "Family", rank_2 = "Class", ...) {
  taxatab <- dplyr::select(
    as.data.frame(unclass(physeq@tax_table)),
    one_of(c(rank_1, rank_2))
  )
  tbl_sum <- gtsummary::tbl_cross(taxatab, ...)
  return(tbl_sum)
}
################################################################################

################################################################################
#' Sankey diagram to compare two taxonomic ranks
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Creates a Sankey (alluvial) diagram to visualize the correspondence between
#' two taxonomic ranks. Useful for comparing taxonomy assignments from different
#' databases or algorithms.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   obtained using the `phyloseq` package.
#' @param rank_1 (character or integer, required) Define the first taxonomic
#'   rank as the number or the name of the column in tax_table slot.
#' @param rank_2 (character or integer, required) Define the second taxonomic
#'   rank as the number or the name of the column in tax_table slot.
#' @param fill_by (character, default "rank_1") Which rank to use for fill
#'   color. Either "rank_1" or "rank_2".
#'
#' @return A ggplot2 object that can be further customized.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' tc_sankey(
#'   Glom_otu,
#'   "Class__eukaryome_Glomero", 
#' "Class"
#' )
#'
#' # Fill by rank_2 instead
#' tc_sankey(
#'   Glom_otu,
#'   "Class__eukaryome_Glomero", 
#'   "Class",
#'   fill_by = "rank_2"
#' )
#'
#' \dontrun{
#' # Add labels to the strata
#' tc_sankey(Glom_otu, "Class__eukaryome_Glomero", "Class") +
#'   geom_label(
#'     stat = "stratum",
#'     aes(label = after_stat(stratum)),
#'     na.rm = TRUE
#'   ) +
#'   theme(legend.position = "none")
#' }
tc_sankey <- function(physeq, rank_1, rank_2, fill_by = "rank_1") {
  rlang::check_installed("ggalluvial", reason = "to create Sankey diagrams")
  verify_pq(physeq)

  fill_by <- match.arg(fill_by, c("rank_1", "rank_2"))
  fill_rank <- if (fill_by == "rank_1") rank_1 else rank_2

  rank_1_name <- colnames(physeq@tax_table)[
    if (is.numeric(rank_1)) rank_1 else which(colnames(physeq@tax_table) == rank_1)
  ]
  rank_2_name <- colnames(physeq@tax_table)[
    if (is.numeric(rank_2)) rank_2 else which(colnames(physeq@tax_table) == rank_2)
  ]
  fill_name <- if (fill_by == "rank_1") rank_1_name else rank_2_name

  df_sank <- as.data.frame(unclass(physeq@tax_table)) |>
    dplyr::select(dplyr::all_of(c(rank_1_name, rank_2_name))) |>
    dplyr::count(.data[[rank_1_name]], .data[[rank_2_name]])

  ggplot(
    df_sank,
    aes(
      axis1 = .data[[rank_1_name]],
      axis2 = .data[[rank_2_name]],
      fill = .data[[fill_name]],
      y = n
    )
  ) +
    ggalluvial::geom_alluvium() +
    ggalluvial::geom_stratum(color = "grey") +
    scale_x_discrete(
      limits = c(rank_1_name, rank_2_name),
      expand = c(.05, .05)
    ) +
    theme_void()
}

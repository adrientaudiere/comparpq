################################################################################
#' Compare two taxonomic assignments of a phyloseq in one call
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Convenience wrapper that assembles, for a single pair of `tax_table` columns
#' (typically the same rank assigned by two databases), the congruence metrics
#' and the four main comparison views of comparpq. It returns a named list
#' bundling [tc_congruence_metrics()], a filtered contingency table, and the
#' [tc_bar()], [tc_sankey()], [tc_heatmap()] plots plus a closure that draws the
#' [tc_circle()] chord diagram on demand (the latter uses base graphics rather
#' than ggplot2).
#'
#' @param physeq (required) A [phyloseq::phyloseq-class()] object.
#' @param rank_1 (character, required) `tax_table` column of the first database
#'   (e.g. `"Genus_SILVA"`).
#' @param rank_2 (character, required) `tax_table` column of the second database
#'   (e.g. `"Genus_KSGP"`).
#' @param rank_3 (character) Column used to color [tc_bar()]. Defaults to the
#'   parent rank of `rank_1` in the same database when it can be resolved from
#'   the column naming (`"<Rank>_<Db>"`), otherwise `rank_1`.
#' @param min_n (integer, default 5) Minimum number of taxa for a contingency
#'   cell to be kept (cells with `n > min_n` are returned).
#'
#' @returns A named list with:
#'   \describe{
#'     \item{congruence}{result of [tc_congruence_metrics()].}
#'     \item{contingency}{data.frame of cells with more than `min_n` taxa,
#'       sorted by decreasing count.}
#'     \item{bar}{a [ggplot2::ggplot()] from [tc_bar()].}
#'     \item{sankey}{a [ggplot2::ggplot()] from [tc_sankey()].}
#'     \item{circle}{a function of no argument that draws [tc_circle()].}
#'     \item{heatmap}{a [ggplot2::ggplot()] from [tc_heatmap()].}
#'   }
#' @export
#' @author Adrien Taudière
#'
#' @seealso [count_taxo_congruence()], [plot_congruence_counts()],
#'   [tc_congruence_metrics()]
#'
#' @examples
#' \donttest{
#' pq <- MiscMetabar::subset_taxa_pq(Glom_otu, phyloseq::taxa_sums(Glom_otu) > 5000)
#' res <- compare_taxo_db(pq, "Genus", "Genus__eukaryome_Glomero")
#' res$contingency
#' res$bar
#' }
compare_taxo_db <- function(physeq, rank_1, rank_2, rank_3 = NULL, min_n = 5) {
  verify_pq(physeq)

  if (is.null(rank_3)) {
    parsed <- split_rank_db(rank_1)
    if (!is.null(parsed)) {
      parent_rank <- rank_before(parsed$rank)
      if (!is.na(parent_rank)) {
        candidate <- paste0(parent_rank, "_", parsed$db)
        if (candidate %in% colnames(physeq@tax_table)) {
          rank_3 <- candidate
        }
      }
    }
  }

  congruence <- tc_congruence_metrics(
    physeq,
    ranks_1 = rank_1,
    ranks_2 = rank_2
  )

  tab <- table(
    physeq@tax_table[, rank_1],
    physeq@tax_table[, rank_2],
    useNA = "ifany"
  )
  contingency <- as.data.frame(tab, stringsAsFactors = FALSE)
  colnames(contingency) <- c(rank_1, rank_2, "n")
  contingency <- contingency[contingency$n > min_n, ]
  contingency <- contingency[order(-contingency$n), ]

  color_rank <- rank_3
  if (is.null(color_rank)) {
    color_rank <- rank_1
  }

  bar_plot <- tc_bar(
    physeq,
    rank_1 = rank_1,
    rank_2 = rank_2,
    color_rank = color_rank
  )

  sankey_plot <- tc_sankey(
    physeq,
    rank_1 = rank_1,
    rank_2 = rank_2,
    fill_by = "rank_1"
  )

  heatmap_plot <- tc_heatmap(
    physeq,
    rank_1 = rank_1,
    rank_2 = rank_2,
    log10trans = TRUE
  )

  circle_fn <- function() {
    tc_circle(physeq, rank_1 = rank_1, rank_2 = rank_2)
  }

  list(
    congruence = congruence,
    contingency = contingency,
    bar = bar_plot,
    sankey = sankey_plot,
    circle = circle_fn,
    heatmap = heatmap_plot
  )
}
################################################################################

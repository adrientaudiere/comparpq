################################################################################
#' Build the grid of pairwise database comparisons across phyloseq objects
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Scans the `tax_table` of each phyloseq object for columns named
#' `"<Rank>_<Db>"` (e.g. `"Genus_SILVA"`, `"Genus_KSGP"`) and enumerates, for
#' each object and each rank, every pairwise combination of databases. The
#' resulting long-format grid has one row per comparison and can be iterated
#' over to run [compare_taxo_db()] or [tc_congruence_metrics()] systematically
#' across objects, ranks and database pairs.
#'
#' @param pq_list (required) A named list of [phyloseq::phyloseq-class()]
#'   objects, or a [list_phyloseq()] object. Names are used in the `physeq`
#'   column of the grid; unnamed lists are auto-named `pq1`, `pq2`, ...
#' @param ranks (character, default `c("Order", "Family", "Genus")`) Taxonomic
#'   ranks to consider. Ranks with fewer than two database columns in a given
#'   phyloseq object are skipped for that object.
#'
#' @returns A [tibble::tibble()] with one row per pairwise comparison and the
#'   columns `physeq`, `rank`, `db_1`, `col_1`, `db_2` and `col_2`, where
#'   `col_1`/`col_2` are the `tax_table` column names to pass to
#'   [compare_taxo_db()]. Returns an empty tibble with the same columns when
#'   no rank carries two databases.
#' @export
#' @author Adrien Taudière
#'
#' @seealso [compare_taxo_db()], [count_taxo_congruence()],
#'   [tc_congruence_metrics()]
#'
#' @examples
#' \donttest{
#' pq <- MiscMetabar::subset_taxa_pq(Glom_otu, phyloseq::taxa_sums(Glom_otu) > 5000)
#' mat <- as(pq@tax_table, "matrix")
#' pq@tax_table <- phyloseq::tax_table(cbind(
#'   mat,
#'   Genus_SILVA = mat[, "Genus"],
#'   Genus_KSGP = mat[, "Genus"]
#' ))
#' build_comparison_grid(list(glom = pq))
#' }
build_comparison_grid <- function(
  pq_list,
  ranks = c("Order", "Family", "Genus")
) {
  if (inherits(pq_list, "list_phyloseq")) {
    pq_list <- pq_list@phyloseq_list
  }
  if (!is.list(pq_list) || length(pq_list) == 0) {
    cli::cli_abort(
      "{.arg pq_list} must be a non-empty named list of phyloseq objects \\
      or a {.cls list_phyloseq} object."
    )
  }
  if (is.null(names(pq_list)) || any(names(pq_list) == "")) {
    names(pq_list) <- paste0("pq", seq_along(pq_list))
    cli::cli_inform(
      "Unnamed {.arg pq_list}: objects auto-named {.val {names(pq_list)}}."
    )
  }
  is_pq <- vapply(pq_list, function(x) inherits(x, "phyloseq"), logical(1))
  if (!all(is_pq)) {
    cli::cli_abort(
      "All elements of {.arg pq_list} must be phyloseq objects; \\
      {.val {names(pq_list)[!is_pq]}} {cli::qty(sum(!is_pq))} {?is/are} not."
    )
  }

  empty_grid <- data.frame(
    physeq = character(0),
    rank = character(0),
    db_1 = character(0),
    col_1 = character(0),
    db_2 = character(0),
    col_2 = character(0),
    stringsAsFactors = FALSE
  )

  rows <- list()
  for (pq_name in names(pq_list)) {
    pq <- pq_list[[pq_name]]
    by_rank <- list()
    for (col in colnames(pq@tax_table)) {
      s <- split_rank_db(col)
      if (is.null(s)) {
        next
      }
      if (!(s$rank %in% ranks)) {
        next
      }
      by_rank[[s$rank]] <- rbind(
        by_rank[[s$rank]],
        data.frame(db = s$db, col = col, stringsAsFactors = FALSE)
      )
    }
    for (rank_name in names(by_rank)) {
      dbs <- by_rank[[rank_name]]
      if (nrow(dbs) < 2) {
        next
      }
      dbs <- dbs[order(dbs$db), ]
      pairs <- utils::combn(seq_len(nrow(dbs)), 2, simplify = FALSE)
      for (p in pairs) {
        rows[[length(rows) + 1]] <- data.frame(
          physeq = pq_name,
          rank = rank_name,
          db_1 = dbs$db[p[1]],
          col_1 = dbs$col[p[1]],
          db_2 = dbs$db[p[2]],
          col_2 = dbs$col[p[2]],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0) {
    return(tibble::as_tibble(empty_grid))
  }
  do.call(rbind, rows) |> tibble::as_tibble()
}
################################################################################

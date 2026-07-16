################################################################################
#' Stacked barplot of taxonomic-assignment congruence across ranks
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Draws, for several taxonomic ranks at once, a stacked barplot of the
#' percentage of taxa (or sequences) falling into each congruence category when
#' comparing two databases (see [count_taxo_congruence()]). The
#' `only_<rank>_<db>` categories are normalized to `only_<db>` so that the color
#' of a database is consistent across ranks (e.g. `only_Order_EUK` and
#' `only_Family_EUK` share the `only_EUK` color).
#'
#' Columns can be supplied in two ways: explicitly through `ranks_1` / `ranks_2`
#' (paired column names), or automatically through `suffix_1` / `suffix_2`
#' (database suffixes such as `"_EUK"`), in which case columns are paired by
#' taxonomic rank.
#'
#' @param physeq (required) A [phyloseq::phyloseq-class()] object.
#' @param ranks_1 (character vector) `tax_table` column names for the first
#'   database, one per comparison. Ignored when `suffix_1` is provided.
#' @param ranks_2 (character vector) `tax_table` column names for the second
#'   database, same length as `ranks_1`. Ignored when `suffix_2` is provided.
#' @param suffix_1 (character) Column suffix (e.g. `"_EUK"`) selecting the
#'   columns of the first database automatically. When provided together with
#'   `suffix_2`, `ranks_1` / `ranks_2` are ignored and columns are paired by
#'   rank.
#' @param suffix_2 (character) Column suffix (e.g. `"_Unite"`) for the second
#'   database.
#' @param ranks (character vector) Rank labels to plot, in the desired y-axis
#'   order (broadest to finest). Accepts standard ranks (`"Species"` ->
#'   `Species_<db>`) and custom prefixes (`"genusSpeciesEpithet"`). A pair is
#'   kept only when both `<prefix><suffix_1>` and `<prefix><suffix_2>` exist;
#'   others are skipped with a warning. `NULL` (default) auto-selects all
#'   standard ranks present in both databases, in taxonomic order.
#' @param n_seq (logical, default `FALSE`) `FALSE`: bars show the percentage of
#'   taxa (`pct_asv`); `TRUE`: bars show the percentage of sequences (`pct_seq`).
#'
#' @returns A [ggplot2::ggplot()] object (stacked bars, one bar per rank, filled
#'   by normalized congruence category).
#' @export
#' @author Adrien Taudière
#'
#' @seealso [count_taxo_congruence()], [compare_taxo_db()]
#'
#' @examples
#' \donttest{
#' plot_congruence_counts(
#'   Glom_otu,
#'   ranks_1 = c("Order", "Family", "Genus"),
#'   ranks_2 = c(
#'     "Order__eukaryome_Glomero",
#'     "Family__eukaryome_Glomero",
#'     "Genus__eukaryome_Glomero"
#'   )
#' )
#' }
plot_congruence_counts <- function(
  physeq,
  ranks_1 = NULL,
  ranks_2 = NULL,
  suffix_1 = NULL,
  suffix_2 = NULL,
  ranks = NULL,
  n_seq = FALSE
) {
  verify_pq(physeq)

  use_suffix <- !is.null(suffix_1) || !is.null(suffix_2)
  rank_labels <- NULL

  if (use_suffix) {
    if (is.null(suffix_1) || is.null(suffix_2)) {
      cli::cli_abort(
        "{.arg suffix_1} and {.arg suffix_2} must be provided together."
      )
    }
    tt_cols <- colnames(physeq@tax_table)

    if (is.null(ranks)) {
      cols_1 <- tt_cols[endsWith(tt_cols, suffix_1)]
      cols_2 <- tt_cols[endsWith(tt_cols, suffix_2)]
      cols_1 <- cols_1[vapply(
        cols_1,
        function(c) !is.null(split_rank_db(c)),
        logical(1)
      )]
      cols_2 <- cols_2[vapply(
        cols_2,
        function(c) !is.null(split_rank_db(c)),
        logical(1)
      )]

      rank_to_col1 <- stats::setNames(
        vapply(cols_1, function(c) split_rank_db(c)$rank, character(1)),
        cols_1
      )
      ordered_ranks <- intersect(known_ranks, rank_to_col1)

      ranks_1 <- ranks_2 <- rank_labels <- character(0)
      for (r in ordered_ranks) {
        c1 <- names(rank_to_col1)[rank_to_col1 == r][1]
        c2_match <- cols_2[vapply(
          cols_2,
          function(c2) {
            s <- split_rank_db(c2)
            !is.null(s) && s$rank == r
          },
          logical(1)
        )]
        if (length(c2_match) > 0) {
          ranks_1 <- c(ranks_1, c1)
          ranks_2 <- c(ranks_2, c2_match[1])
          rank_labels <- c(rank_labels, r)
        }
      }
    } else {
      ranks_1 <- ranks_2 <- rank_labels <- character(0)
      for (rk in ranks) {
        c1 <- paste0(rk, suffix_1)
        c2 <- paste0(rk, suffix_2)
        if (c1 %in% tt_cols && c2 %in% tt_cols) {
          ranks_1 <- c(ranks_1, c1)
          ranks_2 <- c(ranks_2, c2)
          rank_labels <- c(rank_labels, rk)
        } else {
          cli::cli_warn(
            "Rank {.val {rk}} skipped: missing column(s) ({c1} / {c2})."
          )
        }
      }
    }
  } else {
    rank_labels <- vapply(
      ranks_1,
      function(col) {
        parsed <- split_rank_db(col)
        if (!is.null(parsed)) {
          parsed$rank
        } else {
          col
        }
      },
      character(1)
    )
  }

  if (length(ranks_1) != length(ranks_2)) {
    cli::cli_abort(
      "{.arg ranks_1} ({length(ranks_1)}) and {.arg ranks_2} ({length(ranks_2)}) must have the same length."
    )
  }
  if (length(ranks_1) == 0) {
    cli::cli_abort("No column pair to compare. Check {.arg ranks} / suffixes.")
  }

  normalize_cat <- function(cat) {
    if (cat %in% c("both_equal", "both_na", "different")) {
      return(cat)
    }
    if (startsWith(cat, "only_")) {
      col_name <- substr(cat, 6, nchar(cat))
      parsed <- split_rank_db(col_name)
      if (!is.null(parsed)) {
        return(paste0("only_", parsed$db))
      }
      if (!is.null(suffix_1) && endsWith(col_name, suffix_1)) {
        return(paste0("only_", gsub("^_", "", suffix_1)))
      }
      if (!is.null(suffix_2) && endsWith(col_name, suffix_2)) {
        return(paste0("only_", gsub("^_", "", suffix_2)))
      }
    }
    cat
  }

  results <- list()
  for (i in seq_along(ranks_1)) {
    df <- count_taxo_congruence(physeq, ranks_1[i], ranks_2[i])
    df$rank <- rank_labels[i]
    results[[i]] <- df
  }
  all_df <- do.call(rbind, results)
  all_df$category_color <- vapply(all_df$category, normalize_cat, character(1))

  if (use_suffix && !is.null(ranks)) {
    rank_order <- rank_labels
  } else {
    tt_ranks <- vapply(
      colnames(physeq@tax_table),
      function(cn) {
        s <- split_rank_db(cn)
        if (is.null(s)) {
          ""
        } else {
          s$rank
        }
      },
      character(1)
    )
    rank_order <- unique(tt_ranks[tt_ranks != ""])
    extra_ranks <- setdiff(unique(all_df$rank), rank_order)
    if (length(extra_ranks) > 0) {
      sp_idx <- match("Species", rank_order)
      if (!is.na(sp_idx)) {
        rank_order <- append(rank_order, extra_ranks, after = sp_idx)
      } else {
        rank_order <- c(rank_order, extra_ranks)
      }
    }
  }

  all_df$rank <- factor(all_df$rank, levels = rev(rank_order))

  fixed_colors <- c(
    both_equal = "#2e7d32",
    both_na = "#9e9e9e",
    different = "#b3220e"
  )
  only_cats <- unique(all_df$category_color[
    startsWith(all_df$category_color, "only_")
  ])
  only_palette <- c("#1565c0", "#7b1fa2", "#00838f", "#ad1457", "#c46e29")
  only_colors <- stats::setNames(only_palette[seq_along(only_cats)], only_cats)
  color_values <- c(fixed_colors, only_colors)

  cat_order <- c(
    "both_equal",
    "both_na",
    only_cats[order(only_cats)],
    "different"
  )
  all_df$category_color <- factor(all_df$category_color, levels = cat_order)

  y_col <- if (n_seq) {
    "pct_seq"
  } else {
    "pct_asv"
  }
  y_lab <- if (n_seq) {
    "% of sequences"
  } else {
    "% of ASVs"
  }
  all_df$y_value <- all_df[[y_col]]

  ggplot2::ggplot(
    all_df,
    ggplot2::aes(x = .data$rank, y = .data$y_value, fill = .data$category_color)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::scale_fill_manual(values = color_values) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = y_lab, fill = "Category") +
    ggplot2::theme_minimal()
}
################################################################################

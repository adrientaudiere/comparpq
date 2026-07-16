################################################################################
# Internal helpers shared by the taxonomic-database comparison functions
# (`count_taxo_congruence()`, `plot_congruence_counts()`, `compare_taxo_db()`).
# These parse `tax_table` column names of the form "<Rank>_<Db>" and know the
# standard taxonomic hierarchies. They are not exported.
################################################################################

# Recognized taxonomic ranks, sorted by decreasing length so that compound
# names (e.g. "Division_Subdivision") are tested before short prefixes.
known_ranks <- c(
  "Division_Subdivision",
  "Supergroup",
  "Domain",
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

#' Parent rank of a taxonomic rank
#'
#' Returns the rank immediately above `rank_name` in the supported hierarchies
#' (standard Linnaean and PR2), or `NA_character_` when none applies.
#'
#' @param rank_name (character) A taxonomic rank name.
#' @return The parent rank name, or `NA_character_`.
#' @keywords internal
#' @noRd
rank_before <- function(rank_name) {
  hierarchies <- list(
    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    c(
      "Domain",
      "Supergroup",
      "Division_Subdivision",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    )
  )
  for (h in hierarchies) {
    idx <- match(rank_name, h)
    if (!is.na(idx) && idx > 1) {
      return(h[idx - 1])
    }
  }
  NA_character_
}

#' Split a "<Rank>_<Db>" tax_table column name
#'
#' For a column of the form `"<Rank>_<Db>"` returns `list(rank, db)`, or `NULL`
#' when the column is not a rank-by-database column (e.g. blast/RDP helper
#' columns). Correctly handles database names containing an underscore
#' (e.g. `KSGP_archaea`) and compound ranks (e.g. `Division_Subdivision`).
#'
#' @param cn (character) A single column name.
#' @return `list(rank, db)` or `NULL`.
#' @keywords internal
#' @noRd
split_rank_db <- function(cn) {
  if (grepl("^blast_", cn)) {
    return(NULL)
  }
  if (grepl("_RDP_bootstrap$", cn)) {
    return(NULL)
  }
  if (grepl("_RDP$", cn)) {
    return(NULL)
  }
  if (grepl("_BLAST$", cn)) {
    return(NULL)
  }
  for (r in known_ranks) {
    if (startsWith(cn, paste0(r, "_"))) {
      db <- substr(cn, nchar(r) + 2, nchar(cn))
      return(list(rank = r, db = db))
    }
  }
  NULL
}
################################################################################

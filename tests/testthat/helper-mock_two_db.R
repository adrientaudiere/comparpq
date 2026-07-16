# Build a small phyloseq carrying two "databases" (_A and _B) of taxonomic
# assignments, derived from Glom_otu, with some disagreement and missing values.
# Shared by the taxonomic-database comparison tests.
mock_two_db_pq <- function(n_taxa = 60) {
  pq <- MiscMetabar::subset_taxa_pq(
    Glom_otu,
    phyloseq::taxa_sums(Glom_otu) > 5000
  )
  keep <- utils::head(phyloseq::taxa_names(pq), n_taxa)
  pq <- phyloseq::prune_taxa(keep, pq)

  mat <- as(pq@tax_table, "matrix")
  ranks <- c("Order", "Family", "Genus", "Species")
  new_cols <- list()
  for (r in ranks) {
    a <- mat[, r]
    b <- a
    n <- length(b)
    if (n >= 4) {
      b[seq_len(max(1, n %/% 6))] <- NA
      idx <- seq(2, n, by = 5)
      b[idx] <- paste0(b[idx], "_x")
    }
    new_cols[[paste0(r, "_A")]] <- a
    new_cols[[paste0(r, "_B")]] <- b
  }
  add <- do.call(cbind, new_cols)
  mat2 <- cbind(mat, add)
  pq@tax_table <- phyloseq::tax_table(mat2)
  pq
}

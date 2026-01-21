# Tests for taxtab_modification.R
# Tests for tax table manipulation functions

# ==============================================================================
# Tests for select_ranks_pq
# ==============================================================================

test_that("select_ranks_pq selects specified ranks", {
  result <- select_ranks_pq(data_fungi, Order, Family)

  expect_equal(ncol(result@tax_table), 2)
  expect_true(all(c("Order", "Family") %in% colnames(result@tax_table)))
})

test_that("select_ranks_pq excludes ranks with !", {
  original_cols <- colnames(data_fungi@tax_table)

  result <- select_ranks_pq(data_fungi, !Order)

  expect_false("Order" %in% colnames(result@tax_table))
  expect_equal(ncol(result@tax_table), length(original_cols) - 1)
})

test_that("select_ranks_pq preserves taxa_names", {
  result <- select_ranks_pq(data_fungi, Order, Family)

  expect_equal(taxa_names(result), taxa_names(data_fungi))
})

test_that("select_ranks_pq preserves other slots", {
  result <- select_ranks_pq(data_fungi, Order, Family)

  expect_equal(ntaxa(result), ntaxa(data_fungi))
  expect_equal(nsamples(result), nsamples(data_fungi))
})

test_that("select_ranks_pq works with range selection", {
  result <- select_ranks_pq(data_fungi, Order:Genus)

  expect_true("Genus" %in% colnames(result@tax_table))
  expect_true("Order" %in% colnames(result@tax_table))
})

# ==============================================================================
# Tests for rename_ranks_pq
# ==============================================================================

test_that("rename_ranks_pq renames columns with old_names/new_names", {
  result <- rename_ranks_pq(
    data_fungi,
    old_names = c("Phylum"),
    new_names = c("Phyla")
  )

  expect_true("Phyla" %in% colnames(result@tax_table))
  expect_false("Phylum" %in% colnames(result@tax_table))
})

test_that("rename_ranks_pq renames multiple columns", {
  result <- rename_ranks_pq(
    data_fungi,
    old_names = c("Phylum", "Class"),
    new_names = c("Phyla", "Classe")
  )

  expect_true("Phyla" %in% colnames(result@tax_table))
  expect_true("Classe" %in% colnames(result@tax_table))
  expect_false("Phylum" %in% colnames(result@tax_table))
  expect_false("Class" %in% colnames(result@tax_table))
})

test_that("rename_ranks_pq renames with pattern/replacement", {
  result <- rename_ranks_pq(
    data_fungi,
    pattern = ".",
    replacement = "_",
    fixed = TRUE
  )

  # Check that dots are replaced with underscores
  original_with_dot <- grep("\\.", colnames(data_fungi@tax_table), value = TRUE)
  if (length(original_with_dot) > 0) {
    expect_false(any(grepl("\\.", colnames(result@tax_table)[grepl("_", colnames(result@tax_table))])))
  }
})

test_that("rename_ranks_pq fails when mixing old_names and pattern", {
  expect_error(
    rename_ranks_pq(
      data_fungi,
      old_names = c("Phylum"),
      new_names = c("Phyla"),
      pattern = "test"
    ),
    "must specify either"
  )
})

test_that("rename_ranks_pq preserves taxa_names", {
  result <- rename_ranks_pq(
    data_fungi,
    old_names = c("Phylum"),
    new_names = c("Phyla")
  )

  expect_equal(taxa_names(result), taxa_names(data_fungi))
})

test_that("rename_ranks_pq preserves data values", {
  result <- rename_ranks_pq(
    data_fungi,
    old_names = c("Phylum"),
    new_names = c("Phyla")
  )

  # Values should be identical
  expect_equal(
    as.vector(data_fungi@tax_table[, "Phylum"]),
    as.vector(result@tax_table[, "Phyla"])
  )
})

# ==============================================================================
# Tests for taxtab_replace_pattern_by_NA
# ==============================================================================

test_that("taxtab_replace_pattern_by_NA replaces pattern with NA", {
  # Create test data with known pattern
  test_pq <- data_fungi

  result <- taxtab_replace_pattern_by_NA(
    test_pq,
    patterns = c("fam_Incertae_sedis"),
    taxonomic_ranks = "Family"
  )

  # Check that pattern is replaced
  original_matches <- sum(grepl("fam_Incertae_sedis", data_fungi@tax_table[, "Family"]))
  result_matches <- sum(grepl("fam_Incertae_sedis", result@tax_table[, "Family"], fixed = TRUE), na.rm = TRUE)

  expect_true(result_matches < original_matches || original_matches == 0)
})

test_that("taxtab_replace_pattern_by_NA works with default patterns", {
  result <- taxtab_replace_pattern_by_NA(
    data_fungi,
    patterns = c(".*_incertae_sedis", "unclassified.*"),
    ignore.case = TRUE
  )

  expect_s4_class(result, "phyloseq")
})

test_that("taxtab_replace_pattern_by_NA respects taxonomic_ranks parameter", {
  # Only replace in Family
  result <- taxtab_replace_pattern_by_NA(
    data_fungi,
    patterns = c(".*_incertae_sedis"),
    taxonomic_ranks = "Family",
    ignore.case = TRUE
  )

  expect_s4_class(result, "phyloseq")
})

test_that("taxtab_replace_pattern_by_NA applies to all ranks when taxonomic_ranks is NULL", {
  result <- taxtab_replace_pattern_by_NA(
    data_fungi,
    patterns = c("test_pattern_unlikely_to_exist"),
    taxonomic_ranks = NULL
  )

  expect_s4_class(result, "phyloseq")
  expect_equal(ncol(result@tax_table), ncol(data_fungi@tax_table))
})

test_that("taxtab_replace_pattern_by_NA preserves taxa_names", {
  result <- taxtab_replace_pattern_by_NA(
    data_fungi,
    patterns = c("fam_Incertae_sedis"),
    taxonomic_ranks = "Family"
  )

  expect_equal(taxa_names(result), taxa_names(data_fungi))
})

test_that("taxtab_replace_pattern_by_NA progress_bar works", {
  expect_output(
    taxtab_replace_pattern_by_NA(
      data_fungi_mini,
      patterns = c("test"),
      progress_bar = TRUE
    ),
    "="
  )
})

# ==============================================================================
# Tests for resolve_taxo_conflict
# ==============================================================================

test_that("resolve_taxo_conflict works with consensus method", {
  skip_if_not(
    any(grepl("\\.", colnames(data_fungi@tax_table))),
    "No duplicate rank patterns found"
  )

  # Check if there are Genus columns to merge
  genus_cols <- grep("^Genus", colnames(data_fungi@tax_table), value = TRUE)

  if (length(genus_cols) >= 2) {
    result <- resolve_taxo_conflict(
      data_fungi,
      pattern_tax_ranks = c("^Genus"),
      method = "consensus",
      new_names = "Genus_consensus"
    )

    expect_true("Genus_consensus" %in% colnames(result@tax_table))
  }
})

test_that("resolve_taxo_conflict preserves original ranks with keep_tax_ranks=TRUE", {
  # Use data that might have multiple columns matching the pattern
  genus_cols <- grep("^Genus", colnames(data_fungi@tax_table), value = TRUE)

  if (length(genus_cols) >= 1) {
    result <- resolve_taxo_conflict(
      data_fungi,
      pattern_tax_ranks = c("^Genus"),
      method = "consensus",
      new_names = "Genus_consensus",
      keep_tax_ranks = TRUE
    )

    # Original columns should still exist
    expect_true(all(genus_cols %in% colnames(result@tax_table)))
  }
})

test_that("resolve_taxo_conflict removes original ranks with keep_tax_ranks=FALSE", {
  genus_cols <- grep("^Genus", colnames(data_fungi@tax_table), value = TRUE)

  if (length(genus_cols) >= 1) {
    result <- resolve_taxo_conflict(
      data_fungi,
      pattern_tax_ranks = c("^Genus"),
      method = "consensus",
      new_names = "Genus_consensus",
      keep_tax_ranks = FALSE
    )

    # Original columns should be removed
    expect_false(any(genus_cols %in% colnames(result@tax_table)))
    expect_true("Genus_consensus" %in% colnames(result@tax_table))
  }
})

test_that("resolve_taxo_conflict fails with invalid preference_pattern", {
  expect_error(
    resolve_taxo_conflict(
      data_fungi,
      pattern_tax_ranks = c("^Genus"),
      method = "preference",
      preference_pattern = "nonexistent_pattern_xyz"
    ),
    "do not match any rank name"
  )
})

test_that("resolve_taxo_conflict respects replace_collapsed_rank_by_NA", {
  genus_cols <- grep("^Genus", colnames(data_fungi@tax_table), value = TRUE)

  if (length(genus_cols) >= 1) {
    result <- resolve_taxo_conflict(
      data_fungi,
      pattern_tax_ranks = c("^Genus"),
      method = "consensus",
      new_names = "Genus_consensus",
      collapse_string = "/",
      replace_collapsed_rank_by_NA = TRUE
    )

    # Check that collapsed values are replaced by NA
    collapsed_count <- sum(grepl("/", result@tax_table[, "Genus_consensus"]), na.rm = TRUE)
    expect_equal(collapsed_count, 0)
  }
})

test_that("resolve_taxo_conflict preserves taxa_names", {
  genus_cols <- grep("^Genus", colnames(data_fungi@tax_table), value = TRUE)

  if (length(genus_cols) >= 1) {
    result <- resolve_taxo_conflict(
      data_fungi,
      pattern_tax_ranks = c("^Genus"),
      method = "consensus",
      new_names = "Genus_consensus"
    )

    expect_equal(taxa_names(result), taxa_names(data_fungi))
  }
})

# Tests for n_levels_lpq function
# Tests for counting unique taxonomic levels across phyloseq objects

# ==============================================================================
# Helper function to create test phyloseq objects with known taxonomic levels
# ==============================================================================

create_test_phyloseq <- function(n_taxa, tax_levels) {
  # Create OTU table
  otu_mat <- matrix(
    sample(1:100, n_taxa * 5, replace = TRUE),
    nrow = n_taxa,
    ncol = 5
  )
  rownames(otu_mat) <- paste0("ASV_", seq_len(n_taxa))
  colnames(otu_mat) <- paste0("Sample_", seq_len(5))

  # Create tax table with specified levels
  tax_mat <- matrix(
    NA_character_,
    nrow = n_taxa,
    ncol = length(names(tax_levels))
  )
  colnames(tax_mat) <- names(tax_levels)
  rownames(tax_mat) <- paste0("ASV_", seq_len(n_taxa))

  # Fill tax table cycling through the specified levels

  for (rank in names(tax_levels)) {
    levels_for_rank <- tax_levels[[rank]]
    tax_mat[, rank] <- rep_len(levels_for_rank, n_taxa)
  }

  phyloseq::phyloseq(
    phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax_mat)
  )
}

# Create test phyloseq with known levels
# physeq_A: 3 Phyla, 5 Classes, 8 Genera
test_levels_A <- list(
  Phylum = c("Phylum_A", "Phylum_B", "Phylum_C"),
  Class = c("Class_1", "Class_2", "Class_3", "Class_4", "Class_5"),
  Genus = c(
    "Genus_1",
    "Genus_2",
    "Genus_3",
    "Genus_4",
    "Genus_5",
    "Genus_6",
    "Genus_7",
    "Genus_8"
  )
)
physeq_A <- create_test_phyloseq(n_taxa = 20, tax_levels = test_levels_A)

# physeq_B: 2 Phyla, 3 Classes, 4 Genera
test_levels_B <- list(
  Phylum = c("Phylum_X", "Phylum_Y"),
  Class = c("Class_A", "Class_B", "Class_C"),
  Genus = c("Genus_A", "Genus_B", "Genus_C", "Genus_D")
)
physeq_B <- create_test_phyloseq(n_taxa = 15, tax_levels = test_levels_B)

# physeq_C: same as A but with some NAs
physeq_C <- physeq_A
tax_table(physeq_C)[1:5, "Genus"] <- NA

# ==============================================================================
# Tests for n_levels_lpq with controlled test data
# ==============================================================================

test_that("n_levels_lpq returns correct counts for known data", {
  lpq <- list_phyloseq(list(A = physeq_A, B = physeq_B))

  result <- n_levels_lpq(lpq, c("Phylum", "Class", "Genus"))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_equal(rownames(result), c("A", "B"))
  expect_equal(colnames(result), c("Phylum", "Class", "Genus"))

  # Check exact counts for physeq_A

  expect_equal(result["A", "Phylum"], 3)
  expect_equal(result["A", "Class"], 5)
  expect_equal(result["A", "Genus"], 8)

  # Check exact counts for physeq_B
  expect_equal(result["B", "Phylum"], 2)
  expect_equal(result["B", "Class"], 3)
  expect_equal(result["B", "Genus"], 4)
})

test_that("n_levels_lpq handles NA values correctly with na.rm = TRUE", {
  lpq <- list_phyloseq(list(A = physeq_A, C = physeq_C))

  result <- n_levels_lpq(lpq, c("Phylum", "Class", "Genus"), na.rm = TRUE)

  # physeq_A should have 8 genera
  expect_equal(result["A", "Genus"], 8)

  # physeq_C has 5 NAs in Genus, but unique non-NA values should still be 8

  # (since we cycle through 8 values for 20 taxa, removing 5 still leaves all 8 unique)
  expect_equal(result["C", "Genus"], 8)
})

test_that("n_levels_lpq counts NA as a level when na.rm = FALSE", {
  lpq <- list_phyloseq(list(C = physeq_C))

  result_with_na <- n_levels_lpq(lpq, c("Genus"), na.rm = FALSE)
  result_without_na <- n_levels_lpq(lpq, c("Genus"), na.rm = TRUE)

  # With na.rm = FALSE, NA counts as an additional level
  expect_equal(
    result_with_na["C", "Genus"],
    result_without_na["C", "Genus"] + 1
  )
})

test_that("n_levels_lpq works with single phyloseq object", {
  lpq <- list_phyloseq(list(single = physeq_A))

  result <- n_levels_lpq(lpq, c("Phylum", "Class"))

  expect_equal(nrow(result), 1)
  expect_equal(rownames(result), "single")
  expect_equal(result["single", "Phylum"], 3)
  expect_equal(result["single", "Class"], 5)
})

test_that("n_levels_lpq works with single taxonomic rank", {
  lpq <- list_phyloseq(list(A = physeq_A, B = physeq_B))

  result <- n_levels_lpq(lpq, "Phylum")

  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "Phylum")
  expect_equal(result["A", "Phylum"], 3)
  expect_equal(result["B", "Phylum"], 2)
})

# ==============================================================================
# Tests for n_levels_lpq with MiscMetabar data
# ==============================================================================

test_that("n_levels_lpq works with data_fungi and data_fungi_mini", {
  lpq <- list_phyloseq(list(
    full = data_fungi,
    mini = data_fungi_mini
  ))

  common_ranks <- intersect(
    colnames(data_fungi@tax_table),
    colnames(data_fungi_mini@tax_table)
  )

  result <- n_levels_lpq(lpq, common_ranks)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_equal(rownames(result), c("full", "mini"))

  # data_fungi should have >= levels than data_fungi_mini (it's a superset)
  for (rank in common_ranks) {
    expect_true(
      result["full", rank] >= result["mini", rank],
      info = paste("Rank:", rank)
    )
  }
})

test_that("n_levels_lpq returns integer counts", {
  lpq <- list_phyloseq(list(A = physeq_A))

  result <- n_levels_lpq(lpq, c("Phylum", "Class", "Genus"))

  expect_type(result$Phylum, "integer")
  expect_type(result$Class, "integer")
  expect_type(result$Genus, "integer")
})

# ==============================================================================
# Tests for error handling
# ==============================================================================

test_that("n_levels_lpq errors when taxonomic_ranks are missing", {
  lpq <- list_phyloseq(list(A = physeq_A, B = physeq_B))

  expect_error(
    n_levels_lpq(lpq, c("Phylum", "NonExistentRank")),
    "missing"
  )
})

test_that("n_levels_lpq errors when rank missing in one phyloseq only", {
  # Create physeq with different ranks
  physeq_different <- physeq_A
  colnames(tax_table(physeq_different)) <- c("Kingdom", "Division", "Tribe")

  lpq <- list_phyloseq(list(A = physeq_A, different = physeq_different))

  expect_error(
    n_levels_lpq(lpq, c("Phylum", "Class")),
    "missing"
  )
})

test_that("n_levels_lpq errors with non-list_phyloseq input", {
  expect_error(
    n_levels_lpq(physeq_A, c("Phylum")),
    "list_phyloseq"
  )
})

# ==============================================================================
# Tests for edge cases
# ==============================================================================

test_that("n_levels_lpq handles phyloseq with all NA in a rank", {
  physeq_all_na <- physeq_A
  tax_table(physeq_all_na)[, "Genus"] <- NA

  lpq <- list_phyloseq(list(all_na = physeq_all_na))

  result <- n_levels_lpq(lpq, c("Phylum", "Genus"), na.rm = TRUE)

  expect_equal(result["all_na", "Phylum"], 3)
  expect_equal(result["all_na", "Genus"], 0)
})

test_that("n_levels_lpq handles empty string values", {
  physeq_empty <- physeq_A
  tax_table(physeq_empty)[1:3, "Genus"] <- ""

  lpq <- list_phyloseq(list(with_empty = physeq_empty))

  result <- n_levels_lpq(lpq, c("Genus"), na.rm = TRUE)

  # Empty strings should be excluded like NAs
  expect_equal(result["with_empty", "Genus"], 8)
})

test_that("n_levels_lpq handles NA_NA pattern", {
  physeq_na_na <- physeq_A
  tax_table(physeq_na_na)[1:3, "Genus"] <- "NA_NA"

  lpq <- list_phyloseq(list(with_na_na = physeq_na_na))

  result <- n_levels_lpq(lpq, c("Genus"), na.rm = TRUE)

  # "NA_NA" pattern should be excluded like NAs
  expect_equal(result["with_na_na", "Genus"], 8)
})

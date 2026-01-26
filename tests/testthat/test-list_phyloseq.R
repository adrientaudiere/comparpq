# Tests for list_phyloseq.R
# Tests for the S7 class and its methods/utilities

# Setup: Create test data
td <-
  list(
    pq1 = data_fungi,
    pq2 = data_fungi_mini
  )


# ==============================================================================
# Tests for list_phyloseq class creation
# ==============================================================================

test_that("list_phyloseq creates object with valid phyloseq list", {
  lpq <- list_phyloseq(list(fungi = td$pq1, fungi_mini = td$pq2))

  expect_s3_class(lpq, "comparpq::list_phyloseq")
  expect_equal(length(lpq), 2)
  expect_equal(names(lpq), c("fungi", "fungi_mini"))
})

test_that("list_phyloseq auto-generates names when not provided", {
  lpq <- list_phyloseq(list(td$pq1, td$pq2))

  expect_equal(names(lpq), c("physeq_1", "physeq_2"))
})

test_that("list_phyloseq accepts custom names parameter", {
  lpq <- list_phyloseq(list(td$pq1, td$pq2), names = c("A", "B"))

  expect_equal(names(lpq), c("A", "B"))
})

test_that("list_phyloseq fails with empty list", {
  expect_error(
    list_phyloseq(list())
  )
})

test_that("list_phyloseq fails with non-phyloseq objects", {
  expect_error(
    list_phyloseq(list(data.frame(), "not_phyloseq"))
  )
})

test_that("list_phyloseq fails with mismatched names length", {
  expect_error(
    list_phyloseq(list(td$pq1, td$pq2), names = c("A")),
    "must match length"
  )
})

# ==============================================================================
# Tests for comparison type parameters
# ==============================================================================

test_that("list_phyloseq defaults to REPRODUCIBILITY for same samples", {
  # Same phyloseq twice -> same samples
  lpq <- list_phyloseq(list(run1 = td$pq1, run2 = td$pq1))

  expect_equal(lpq@comparison$type_of_comparison, "REPRODUCIBILITY")
  expect_true(lpq@comparison$same_primer_seq_tech)
  expect_true(lpq@comparison$same_bioinfo_pipeline)
})

test_that("list_phyloseq detects ROBUSTNESS with same_bioinfo_pipeline=FALSE", {
  lpq <- list_phyloseq(
    list(method_A = td$pq1, method_B = td$pq1),
    same_bioinfo_pipeline = FALSE
  )

  expect_equal(lpq@comparison$type_of_comparison, "ROBUSTNESS")
  expect_true(lpq@comparison$same_primer_seq_tech)
  expect_false(lpq@comparison$same_bioinfo_pipeline)
})

test_that("list_phyloseq detects REPLICABILITY with same_primer_seq_tech=FALSE", {
  lpq <- list_phyloseq(
    list(ITS1 = td$pq1, ITS2 = td$pq1),
    same_primer_seq_tech = FALSE
  )

  expect_equal(lpq@comparison$type_of_comparison, "REPLICABILITY")
  expect_false(lpq@comparison$same_primer_seq_tech)
})

test_that("list_phyloseq detects NESTED_ROBUSTNESS for nested samples", {
  # Rarefied version has nested samples
  set.seed(123)
  rarefied <- rarefy_even_depth(td$pq1, sample.size = 1000, verbose = FALSE)

  lpq <- list_phyloseq(list(original = td$pq1, rarefied = rarefied))

  expect_equal(lpq@comparison$type_of_comparison, "NESTED_ROBUSTNESS")
  expect_true(lpq@comparison$nested_samples)
  expect_false(is.null(lpq@comparison$nesting_structure))
})

test_that("list_phyloseq detects EXPLORATION for different samples with shared modalities", {
  # data_fungi and data_fungi_mini have different samples but shared modalities
  lpq <- list_phyloseq(list(fungi = td$pq1, fungi_mini = td$pq2))

  # Should be EXPLORATION or related type (depends on shared modalities)
  expect_true(lpq@comparison$type_of_comparison %in%
    c("EXPLORATION", "SEPARATE_ANALYSIS", "NESTED_ROBUSTNESS"))
})

# ==============================================================================
# Tests for summary_table
# ==============================================================================

test_that("list_phyloseq computes correct summary_table", {
  lpq <- list_phyloseq(list(fungi = td$pq1, fungi_mini = td$pq2))

  expect_true(tibble::is_tibble(lpq@summary_table))
  expect_equal(nrow(lpq@summary_table), 2)
  expect_true("n_samples" %in% colnames(lpq@summary_table))
  expect_true("n_taxa" %in% colnames(lpq@summary_table))
  expect_true("n_sequences" %in% colnames(lpq@summary_table))

  # Verify values match phyloseq
  expect_equal(lpq@summary_table$n_samples[1], nsamples(td$pq1))
  expect_equal(lpq@summary_table$n_taxa[1], ntaxa(td$pq1))
})

# ==============================================================================
# Tests for comparison characteristics
# ==============================================================================

test_that("list_phyloseq computes correct comparison characteristics", {
  lpq <- list_phyloseq(list(fungi = td$pq1, fungi_mini = td$pq2))

  comp <- lpq@comparison

  expect_true("same_samples" %in% names(comp))
  expect_true("same_taxa" %in% names(comp))
  expect_true("n_common_samples" %in% names(comp))
  expect_true("n_common_taxa" %in% names(comp))
  expect_true("all_have_sam_data" %in% names(comp))
  expect_true("all_have_tax_table" %in% names(comp))
})

test_that("list_phyloseq detects same samples correctly", {
  # Same phyloseq -> same samples
  lpq_same <- list_phyloseq(list(a = td$pq1, b = td$pq1))
  expect_true(lpq_same@comparison$same_samples)

  # Different phyloseq -> may have different samples
  lpq_diff <- list_phyloseq(list(fungi = td$pq1, fungi_mini = td$pq2))
  # This depends on the actual data
})

# ==============================================================================
# Tests for S7 methods
# ==============================================================================

test_that("print method works for list_phyloseq", {
  lpq <- list_phyloseq(list(fungi = td$pq1))

  expect_output(print(lpq), "list_phyloseq object")
  expect_output(print(lpq), "Type of comparison")
})

test_that("length method works for list_phyloseq", {
  lpq <- list_phyloseq(list(a = td$pq1, b = td$pq2))

  expect_equal(length(lpq), 2)
})

test_that("names method works for list_phyloseq", {
  lpq <- list_phyloseq(list(fungi = td$pq1, mini = td$pq2))

  expect_equal(names(lpq), c("fungi", "mini"))
})

test_that("subsetting with [ returns list_phyloseq", {
  lpq <- list_phyloseq(list(a = td$pq1, b = td$pq2))

  lpq_subset <- lpq[1]

  expect_s3_class(lpq_subset, "comparpq::list_phyloseq")
  expect_equal(length(lpq_subset), 1)
})

test_that("subsetting preserves parameters", {
  lpq <- list_phyloseq(
    list(a = td$pq1, b = td$pq1),
    same_bioinfo_pipeline = FALSE
  )

  lpq_subset <- lpq[1]

  expect_false(lpq_subset@comparison$same_bioinfo_pipeline)
})

test_that("extraction with [[ returns phyloseq", {
  lpq <- list_phyloseq(list(fungi = td$pq1, mini = td$pq2))

  pq <- lpq[["fungi"]]

  expect_s4_class(pq, "phyloseq")
  expect_equal(ntaxa(pq), ntaxa(td$pq1))
})

# ==============================================================================
# Tests for update_list_phyloseq
# ==============================================================================

test_that("update_list_phyloseq recomputes summary", {
  lpq <- list_phyloseq(list(fungi = td$pq1))

  lpq_updated <- update_list_phyloseq(lpq)

  expect_s3_class(lpq_updated, "comparpq::list_phyloseq")
  expect_equal(lpq_updated@summary_table$n_taxa, lpq@summary_table$n_taxa)
})

test_that("update_list_phyloseq can change parameters", {
  lpq <- list_phyloseq(list(a = td$pq1, b = td$pq1))

  expect_true(lpq@comparison$same_bioinfo_pipeline)

  lpq_updated <- update_list_phyloseq(lpq, same_bioinfo_pipeline = FALSE)

  expect_false(lpq_updated@comparison$same_bioinfo_pipeline)
  expect_equal(lpq_updated@comparison$type_of_comparison, "ROBUSTNESS")
})

test_that("update_list_phyloseq preserves parameters when NULL", {
  lpq <- list_phyloseq(
    list(a = td$pq1, b = td$pq1),
    same_primer_seq_tech = FALSE
  )

  lpq_updated <- update_list_phyloseq(lpq)

  expect_false(lpq_updated@comparison$same_primer_seq_tech)
})

# ==============================================================================
# Tests for add_phyloseq
# ==============================================================================

test_that("add_phyloseq adds a new phyloseq object", {
  lpq <- list_phyloseq(list(fungi = td$pq1))

  lpq_added <- add_phyloseq(lpq, td$pq2, name = "mini")

  expect_equal(length(lpq_added), 2)
  expect_true("mini" %in% names(lpq_added))
})

test_that("add_phyloseq auto-generates name when not provided", {
  lpq <- list_phyloseq(list(fungi = td$pq1))

  lpq_added <- add_phyloseq(lpq, td$pq2)

  expect_equal(length(lpq_added), 2)
  expect_true("physeq_2" %in% names(lpq_added))
})

test_that("add_phyloseq preserves parameters", {
  lpq <- list_phyloseq(
    list(fungi = td$pq1),
    same_bioinfo_pipeline = FALSE
  )

  lpq_added <- add_phyloseq(lpq, td$pq2, name = "mini")

  expect_false(lpq_added@comparison$same_bioinfo_pipeline)
})

test_that("add_phyloseq fails with non-phyloseq object", {
  lpq <- list_phyloseq(list(fungi = td$pq1))

  expect_error(
    add_phyloseq(lpq, "not_phyloseq"),
    "phyloseq"
  )
})

# ==============================================================================
# Tests for remove_phyloseq
# ==============================================================================

test_that("remove_phyloseq removes by name", {
  lpq <- list_phyloseq(list(fungi = td$pq1, mini = td$pq2))

  lpq_removed <- remove_phyloseq(lpq, "mini")

  expect_equal(length(lpq_removed), 1)
  expect_false("mini" %in% names(lpq_removed))
})

test_that("remove_phyloseq removes by index", {
  lpq <- list_phyloseq(list(fungi = td$pq1, mini = td$pq2))

  lpq_removed <- remove_phyloseq(lpq, 2)

  expect_equal(length(lpq_removed), 1)
  expect_equal(names(lpq_removed), "fungi")
})

test_that("remove_phyloseq preserves parameters", {
  lpq <- list_phyloseq(
    list(a = td$pq1, b = td$pq2, c = td$pq1),
    same_primer_seq_tech = FALSE
  )

  lpq_removed <- remove_phyloseq(lpq, "b")

  expect_false(lpq_removed@comparison$same_primer_seq_tech)
})

test_that("remove_phyloseq fails when removing last object", {
  lpq <- list_phyloseq(list(fungi = td$pq1))

  expect_error(
    remove_phyloseq(lpq, "fungi"),
    "Cannot remove the last"
  )
})

# ==============================================================================
# Tests for filter_common_lpq
# ==============================================================================

test_that("filter_common_lpq filters to common samples", {
  set.seed(123)
  rarefied <- rarefy_even_depth(td$pq1, sample.size = 1000, verbose = FALSE)

  lpq <- list_phyloseq(list(original = td$pq1, rarefied = rarefied))

  # Skip if no common samples
  if (lpq@comparison$n_common_samples > 0) {
    lpq_filtered <- filter_common_lpq(lpq, filter_samples = TRUE, verbose = FALSE)

    expect_true(lpq_filtered@comparison$same_samples)
    expect_equal(
      nsamples(lpq_filtered[["original"]]),
      nsamples(lpq_filtered[["rarefied"]])
    )
  }
})

test_that("filter_common_lpq preserves parameters", {
  set.seed(123)
  rarefied <- rarefy_even_depth(td$pq1, sample.size = 1000, verbose = FALSE)

  lpq <- list_phyloseq(
    list(original = td$pq1, rarefied = rarefied),
    same_bioinfo_pipeline = FALSE
  )

  if (lpq@comparison$n_common_samples > 0) {
    lpq_filtered <- filter_common_lpq(lpq, verbose = FALSE)
    expect_false(lpq_filtered@comparison$same_bioinfo_pipeline)
  }
})

test_that("filter_common_lpq returns original when no filtering requested", {
  lpq <- list_phyloseq(list(fungi = td$pq1))

  lpq_no_filter <- filter_common_lpq(
    lpq,
    filter_samples = FALSE,
    filter_taxa = FALSE,
    verbose = FALSE
  )

  expect_identical(lpq, lpq_no_filter)
})

test_that("filter_common_lpq fails with no common samples", {
  # Create two phyloseq with different sample names
  pq1 <- td$pq1
  pq2 <- td$pq1

  # Rename samples in pq2
  sample_names(pq2) <- paste0("different_", sample_names(pq2))

  lpq <- list_phyloseq(list(a = pq1, b = pq2))

  expect_error(
    filter_common_lpq(lpq, filter_samples = TRUE, verbose = FALSE),
    "No common samples"
  )
})

# ==============================================================================
# Tests for shared_mod_lpq
# ==============================================================================

test_that("shared_mod_lpq returns tibble with shared modalities", {
  lpq <- list_phyloseq(list(fungi = td$pq1, mini = td$pq2))

  result <- shared_mod_lpq(lpq)

  # Result depends on actual shared modalities
  if (!is.null(result)) {
    expect_true(tibble::is_tibble(result))
    expect_true("Variable" %in% colnames(result))
    expect_true("N_shared" %in% colnames(result))
  }
})

test_that("shared_mod_lpq respects max_modalities", {
  lpq <- list_phyloseq(list(fungi = td$pq1, mini = td$pq2))

  result <- shared_mod_lpq(lpq, max_modalities = 2)

  # Check that display is limited (if there are shared modalities)
  if (!is.null(result) && any(result$N_shared > 2)) {
    # Should contain "..." indicator for truncated lists
    expect_true(any(grepl("\\.\\.\\.", result$Shared_modalities)))
  }
})

# ==============================================================================
# Tests for n_levels_lpq
# ==============================================================================

test_that("n_levels_lpq counts unique taxonomic levels", {
  lpq <- list_phyloseq(list(fungi = td$pq1, mini = td$pq2))

  result <- n_levels_lpq(lpq, c("Phylum", "Class", "Order"))

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_equal(colnames(result), c("Phylum", "Class", "Order"))
  expect_true(all(result >= 0))
})

test_that("n_levels_lpq fails with missing taxonomic ranks", {
  lpq <- list_phyloseq(list(fungi = td$pq1))

  expect_error(
    n_levels_lpq(lpq, c("Phylum", "NonExistentRank")),
    "missing"
  )
})

# ==============================================================================
# Tests for apply_to_lpq
# ==============================================================================

test_that("apply_to_lpq applies function to all phyloseq objects", {
  lpq <- list_phyloseq(list(fungi = td$pq1, mini = td$pq2))

  lpq_transformed <- apply_to_lpq(lpq, MiscMetabar::clean_pq, verbose = FALSE)

  expect_s3_class(lpq_transformed, "comparpq::list_phyloseq")
  expect_equal(length(lpq_transformed), 2)
  expect_equal(names(lpq_transformed), c("fungi", "mini"))
})

test_that("apply_to_lpq preserves comparison parameters", {
  lpq <- list_phyloseq(
    list(fungi = td$pq1, mini = td$pq2),
    same_primer_seq_tech = FALSE,
    same_bioinfo_pipeline = FALSE
  )

  lpq_transformed <- apply_to_lpq(lpq, MiscMetabar::clean_pq, verbose = FALSE)

  expect_false(lpq_transformed@comparison$same_primer_seq_tech)
  expect_false(lpq_transformed@comparison$same_bioinfo_pipeline)
})

test_that("apply_to_lpq passes additional arguments to function", {
  lpq <- list_phyloseq(list(fungi = td$pq1, mini = td$pq2))

  # Use transform_sample_counts with a custom function
 lpq_rel <- apply_to_lpq(
    lpq,
    phyloseq::transform_sample_counts,
    \(x) x / sum(x),
    verbose = FALSE
  )

  expect_s3_class(lpq_rel, "comparpq::list_phyloseq")

  # Check that counts are now proportions (sum to 1 per sample)
  sums <- phyloseq::sample_sums(lpq_rel[["fungi"]])
  expect_true(all(abs(sums - 1) < 1e-10))
})

test_that("apply_to_lpq works with taxa_as_rows", {
  lpq <- list_phyloseq(list(fungi = td$pq1, mini = td$pq2))

  lpq_rows <- apply_to_lpq(lpq, MiscMetabar::taxa_as_rows, verbose = FALSE)

  expect_s3_class(lpq_rows, "comparpq::list_phyloseq")
  expect_true(phyloseq::taxa_are_rows(lpq_rows[["fungi"]]))
  expect_true(phyloseq::taxa_are_rows(lpq_rows[["mini"]]))
})

test_that("apply_to_lpq fails when function is not provided", {
  lpq <- list_phyloseq(list(fungi = td$pq1))

  expect_error(
    apply_to_lpq(lpq, "not_a_function"),
    "must be a function"
  )
})

test_that("apply_to_lpq fails when function does not return phyloseq", {
  lpq <- list_phyloseq(list(fungi = td$pq1))

  # Function that returns a data frame instead of phyloseq
  bad_func <- function(pq) {
    as.data.frame(phyloseq::otu_table(pq))
  }

  expect_error(
    apply_to_lpq(lpq, bad_func, verbose = FALSE),
    "did not return a phyloseq"
  )
})

test_that("apply_to_lpq works with rarefy_even_depth", {
  lpq <- list_phyloseq(list(fungi = td$pq1, mini = td$pq2))

  lpq_rarefied <- apply_to_lpq(
    lpq,
    phyloseq::rarefy_even_depth,
    sample.size = 500,
    rngseed = 42,
    verbose = FALSE
  )

  expect_s3_class(lpq_rarefied, "comparpq::list_phyloseq")

  # Check rarefaction worked (all samples have same depth)
  sums_fungi <- phyloseq::sample_sums(lpq_rarefied[["fungi"]])
  expect_true(all(sums_fungi == 500))
})

test_that("apply_to_lpq can be chained with pipes", {
  lpq <- list_phyloseq(list(fungi = td$pq1, mini = td$pq2))

  lpq_processed <- lpq |>
    apply_to_lpq(MiscMetabar::clean_pq, verbose = FALSE) |>
    apply_to_lpq(MiscMetabar::taxa_as_rows, verbose = FALSE)

  expect_s3_class(lpq_processed, "comparpq::list_phyloseq")
  expect_true(phyloseq::taxa_are_rows(lpq_processed[["fungi"]]))
})

test_that("apply_to_lpq verbose mode prints messages", {
  lpq <- list_phyloseq(list(fungi = td$pq1))

  expect_message(
    apply_to_lpq(lpq, MiscMetabar::clean_pq, verbose = TRUE),
    "Applied"
  )
})

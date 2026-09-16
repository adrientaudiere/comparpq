# Tests for compare_taxo.R
# Tests for accuracy metrics comparing taxonomic assignments

# ==============================================================================
# Helper function to create test data with fake taxa
# ==============================================================================

setup_mock_data <- function() {
  # Add fake sequences to data_fungi_mini
  set.seed(123)
  pq_with_fake <- add_shuffle_seq_pq(data_fungi_mini, n_fake = 10)

  # Get true values from original taxonomy
  true_genera <- unique(na.omit(as.vector(data_fungi_mini@tax_table[,
    "Genus"
  ])))
  true_families <- unique(na.omit(as.vector(data_fungi_mini@tax_table[,
    "Family"
  ])))

  list(
    physeq = pq_with_fake,
    true_genera = true_genera,
    true_families = true_families
  )
}

# ==============================================================================
# Tests for tc_metrics_mock_vec
# ==============================================================================

test_that("tc_metrics_mock_vec returns list with correct metrics when fake_taxa=TRUE", {
  td <- setup_mock_data()

  result <- tc_metrics_mock_vec(
    td$physeq,
    taxonomic_rank = "Genus",
    true_values = td$true_genera,
    fake_taxa = TRUE,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_true("TP" %in% names(result))
  expect_true("FP" %in% names(result))
  expect_true("FN" %in% names(result))
  expect_true("TN" %in% names(result))
  expect_true("FDR" %in% names(result))
  expect_true("TPR" %in% names(result))
  expect_true("TNR" %in% names(result))
  expect_true("PPV" %in% names(result))
  expect_true("F1_score" %in% names(result))
  expect_true("MCC" %in% names(result))
  expect_true("ACC" %in% names(result))
})

test_that("tc_metrics_mock_vec returns list with basic metrics when fake_taxa=FALSE", {
  result <- tc_metrics_mock_vec(
    data_fungi_mini,
    taxonomic_rank = "Genus",
    true_values = unique(na.omit(as.vector(data_fungi_mini@tax_table[,
      "Genus"
    ]))),
    fake_taxa = FALSE,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_true("TP" %in% names(result))
  expect_true("FP" %in% names(result))
  expect_true("FN" %in% names(result))
  expect_true("FDR" %in% names(result))
  expect_true("TPR" %in% names(result))
  expect_true("PPV" %in% names(result))
  expect_true("F1_score" %in% names(result))

  # TN, MCC, ACC should NOT be present without fake taxa

  expect_false("TN" %in% names(result))
  expect_false("MCC" %in% names(result))
  expect_false("ACC" %in% names(result))
})

test_that("tc_metrics_mock_vec TP/FP/FN are non-negative integers", {
  td <- setup_mock_data()

  result <- tc_metrics_mock_vec(
    td$physeq,
    taxonomic_rank = "Genus",
    true_values = td$true_genera,
    fake_taxa = TRUE,
    verbose = FALSE
  )

  expect_true(result$TP >= 0)
  expect_true(result$FP >= 0)
  expect_true(result$FN >= 0)
  expect_true(result$TN >= 0)
})

test_that("tc_metrics_mock_vec MCC does not overflow with large counts", {
  n <- 300
  genus <- c(
    paste0("true_", seq_len(n)),
    paste0("wrong_", seq_len(n)),
    rep(NA, n)
  )
  taxa <- c(paste0("asv_", seq_len(2 * n)), paste0("fake_", seq_len(n)))
  otu <- matrix(
    1L,
    nrow = 3 * n,
    ncol = 2,
    dimnames = list(taxa, c("s1", "s2"))
  )
  tax <- matrix(
    c(rep("Fungi", 3 * n), genus),
    ncol = 2,
    dimnames = list(taxa, c("Kingdom", "Genus"))
  )
  physeq <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(tax))
  true_values <- c(paste0("true_", seq_len(n)), paste0("missing_", seq_len(n)))

  result <- expect_no_warning(tc_metrics_mock_vec(
    physeq,
    taxonomic_rank = "Genus",
    true_values = true_values,
    verbose = FALSE
  ))

  expect_equal(c(result$TP, result$FP, result$FN, result$TN), c(n, n, n, n))
  expect_equal(result$MCC, 0)
})

test_that("tc_metrics_mock_vec MCC is 0 when a confusion-matrix margin is 0", {
  taxa <- c("asv_1", "asv_2", "fake_1", "fake_2")
  otu <- matrix(1L, nrow = 4, ncol = 1, dimnames = list(taxa, "s1"))
  tax <- matrix(
    c(rep("Fungi", 4), NA, NA, "G1", "G2"),
    ncol = 2,
    dimnames = list(taxa, c("Kingdom", "Genus"))
  )
  physeq <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(tax))

  unassigned <- tc_metrics_mock_vec(
    physeq,
    taxonomic_rank = "Genus",
    true_values = c("G1", "G3"),
    verbose = FALSE
  )
  expect_equal(unassigned$TP + unassigned$FP, 0)
  expect_equal(unassigned$MCC, 0)

  all_controls_assigned <- tc_metrics_mock_vec(
    physeq,
    taxonomic_rank = "Kingdom",
    true_values = "Fungi",
    verbose = FALSE
  )
  expect_equal(all_controls_assigned$TN + all_controls_assigned$FN, 0)
  expect_equal(all_controls_assigned$MCC, 0)
})

test_that("tc_metrics_mock_vec rates are between 0 and 1", {
  td <- setup_mock_data()

  result <- tc_metrics_mock_vec(
    td$physeq,
    taxonomic_rank = "Genus",
    true_values = td$true_genera,
    fake_taxa = TRUE,
    verbose = FALSE
  )

  # FDR can be NaN if FP + TP = 0
  if (!is.nan(result$FDR)) {
    expect_true(result$FDR >= 0 && result$FDR <= 1)
  }

  # TPR can be NaN if TP + FN = 0
  if (!is.nan(result$TPR)) {
    expect_true(result$TPR >= 0 && result$TPR <= 1)
  }

  # PPV can be NaN if TP + FP = 0
  if (!is.nan(result$PPV)) {
    expect_true(result$PPV >= 0 && result$PPV <= 1)
  }

  # ACC should be between 0 and 1
  if (!is.nan(result$ACC)) {
    expect_true(result$ACC >= 0 && result$ACC <= 1)
  }
})

test_that("tc_metrics_mock_vec respects custom fake_pattern", {
  td <- setup_mock_data()

  result <- tc_metrics_mock_vec(
    td$physeq,
    taxonomic_rank = "Genus",
    true_values = td$true_genera,
    fake_taxa = TRUE,
    fake_pattern = c("^fake_"),
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_true(result$TN >= 0)
})

test_that("tc_metrics_mock_vec verbose outputs message", {
  td <- setup_mock_data()

  expect_message(
    tc_metrics_mock_vec(
      td$physeq,
      taxonomic_rank = "Genus",
      true_values = td$true_genera,
      fake_taxa = TRUE,
      verbose = TRUE
    ),
    "fake taxa were found"
  )
})

test_that("tc_metrics_mock_vec works with integer taxonomic_rank", {
  td <- setup_mock_data()

  # Get the column index for Genus
  genus_idx <- which(colnames(td$physeq@tax_table) == "Genus")

  result <- tc_metrics_mock_vec(
    td$physeq,
    taxonomic_rank = genus_idx,
    true_values = td$true_genera,
    fake_taxa = TRUE,
    verbose = FALSE
  )

  expect_type(result, "list")
})

# ==============================================================================
# Tests for tc_metrics_mock
# ==============================================================================

test_that("tc_metrics_mock returns dataframe with correct structure", {
  td <- setup_mock_data()

  # Create ranks_df with one method and two ranks
  ranks_df <- data.frame(
    method1 = c("Genus", "Family")
  )

  # Create true_values_df
  true_values_df <- data.frame(
    Genus = td$true_genera[1:min(5, length(td$true_genera))],
    Family = td$true_families[1:min(5, length(td$true_families))]
  )

  result <- tc_metrics_mock(
    td$physeq,
    ranks_df = ranks_df,
    true_values_df = true_values_df,
    fake_taxa = TRUE,
    fake_pattern = c("^fake_")
  )

  expect_s3_class(result, "data.frame")
  expect_equal(
    colnames(result),
    c("method_db", "tax_level", "metrics", "values")
  )
})

test_that("tc_metrics_mock fails with mismatched dimensions", {
  td <- setup_mock_data()

  # Create ranks_df with 2 ranks
  ranks_df <- data.frame(
    method1 = c("Genus", "Family")
  )

  # Create true_values_df with only 1 column (mismatch)
  true_values_df <- data.frame(
    Genus = td$true_genera[1:5]
  )

  expect_error(
    tc_metrics_mock(
      td$physeq,
      ranks_df = ranks_df,
      true_values_df = true_values_df
    ),
    "number of rows of ranks_df must be equal to the number of column"
  )
})

test_that("tc_metrics_mock works with multiple methods", {
  td <- setup_mock_data()

  # Create ranks_df with two methods
  ranks_df <- data.frame(
    method_A = c("Genus"),
    method_B = c("Genus")
  )

  # Create true_values_df
  true_values_df <- data.frame(
    Genus = td$true_genera[1:min(5, length(td$true_genera))]
  )

  result <- tc_metrics_mock(
    td$physeq,
    ranks_df = ranks_df,
    true_values_df = true_values_df,
    fake_taxa = TRUE,
    fake_pattern = c("^fake_")
  )

  expect_s3_class(result, "data.frame")
  expect_true("method_A" %in% result$method_db)
  expect_true("method_B" %in% result$method_db)
})

test_that("tc_metrics_mock respects fake_taxa parameter", {
  # Without fake taxa
  result_no_fake <- tc_metrics_mock(
    data_fungi_mini,
    ranks_df = data.frame(method1 = "Genus"),
    true_values_df = data.frame(
      Genus = unique(na.omit(as.vector(data_fungi_mini@tax_table[, "Genus"])))[
        1:5
      ]
    ),
    fake_taxa = FALSE
  )

  # TN, MCC, ACC should not be in metrics
  expect_false("TN" %in% result_no_fake$metrics)
  expect_false("MCC" %in% result_no_fake$metrics)
  expect_false("ACC" %in% result_no_fake$metrics)
})

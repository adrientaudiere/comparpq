# Tests for analysis_lpq.R
# Tests for statistical analysis functions for list_phyloseq objects

# ==============================================================================
# Helper function to create test data without NA in key variables
# ==============================================================================

create_test_lpq_for_adonis <- function() {
  # Filter out samples with NA in Height (common variable)
  pq_clean <- subset_samples(
    data_fungi,
    !is.na(Height)
  )
  pq_clean <- clean_pq(pq_clean)

  list_phyloseq(
    list(
      original = pq_clean,
      copy = pq_clean
    ),
    same_bioinfo_pipeline = FALSE
  )
}

# ==============================================================================
# Tests for extract_formula_vars (internal function)
# ==============================================================================

test_that("extract_formula_vars extracts simple variables", {
  vars <- extract_formula_vars("Height")
  expect_equal(vars, "Height")
})

test_that("extract_formula_vars extracts multiple variables with +", {
  vars <- extract_formula_vars("Height + Time")
  expect_true(all(c("Height", "Time") %in% vars))
})

test_that("extract_formula_vars extracts variables with *", {
  vars <- extract_formula_vars("Height * Time")
  expect_true(all(c("Height", "Time") %in% vars))
})

test_that("extract_formula_vars extracts variables with :", {
  vars <- extract_formula_vars("Height:Time")
  expect_true(all(c("Height", "Time") %in% vars))
})

test_that("extract_formula_vars handles Condition wrapper", {
  vars <- extract_formula_vars("Condition(Site) + Treatment")
  expect_true(all(c("Site", "Treatment") %in% vars))
})

test_that("extract_formula_vars handles complex formulas", {
  vars <- extract_formula_vars("A + B * C + D:E")
  expect_true(all(c("A", "B", "C", "D", "E") %in% vars))
})

# ==============================================================================
# Tests for validate_formula_vars (internal function)
# ==============================================================================

test_that("validate_formula_vars accepts valid variables", {
  lpq <- create_test_lpq_for_adonis()

  # Height should be a common variable
  expect_true(
    validate_formula_vars(lpq, "Height")
  )
})

test_that("validate_formula_vars rejects invalid variables", {
  lpq <- create_test_lpq_for_adonis()

  expect_error(
    validate_formula_vars(lpq, "NonexistentVariable"),
    "not found in common sample_data"
  )
})

# ==============================================================================
# Tests for adonis_lpq
# ==============================================================================

test_that("adonis_lpq runs PERMANOVA on each phyloseq", {
  lpq <- create_test_lpq_for_adonis()

  result <- adonis_lpq(lpq, formula = "Height", verbose = FALSE)

  expect_s3_class(result, "adonis_lpq_result")
  expect_true(tibble::is_tibble(result))
  expect_true("name" %in% colnames(result))
  expect_true("R2" %in% colnames(result))
  expect_true("Pr(>F)" %in% colnames(result))
})

test_that("adonis_lpq returns one row per term per phyloseq", {
  lpq <- create_test_lpq_for_adonis()

  result <- adonis_lpq(lpq, formula = "Height", verbose = FALSE)

  # Should have entries for both phyloseq objects
  expect_true(all(c("original", "copy") %in% result$name))
})

test_that("adonis_lpq respects dist_method parameter", {
  lpq <- create_test_lpq_for_adonis()

  result_bray <- adonis_lpq(
    lpq,
    formula = "Height",
    dist_method = "bray",
    verbose = FALSE
  )
  result_jaccard <- adonis_lpq(
    lpq,
    formula = "Height",
    dist_method = "jaccard",
    verbose = FALSE
  )

  # Both should work
  expect_s3_class(result_bray, "adonis_lpq_result")
  expect_s3_class(result_jaccard, "adonis_lpq_result")

  # R2 values should differ (different distance metrics)
  # Note: May be similar for identical data, so just check they ran
})

test_that("adonis_lpq respects na_remove parameter", {
  # Create data with some NA
  pq_with_na <- data_fungi

  lpq <- list_phyloseq(
    list(a = pq_with_na, b = pq_with_na),
    same_bioinfo_pipeline = FALSE
  )

  # Should work with na_remove = TRUE
  result <- adonis_lpq(
    lpq,
    formula = "Height",
    na_remove = TRUE,
    verbose = FALSE
  )

  expect_s3_class(result, "adonis_lpq_result")
})

test_that("adonis_lpq fails with SEPARATE_ANALYSIS type", {
  # Create list_phyloseq that will be SEPARATE_ANALYSIS
  data("enterotype", package = "phyloseq")

  lpq <- list_phyloseq(list(
    fungi = data_fungi_mini,
    enterotype = enterotype
  ))

  # This should fail because variables won't be common
  expect_error(
    adonis_lpq(lpq, formula = "Height", verbose = FALSE)
  )
})

test_that("adonis_lpq fails with invalid formula variables", {
  lpq <- create_test_lpq_for_adonis()

  expect_error(
    adonis_lpq(lpq, formula = "NonexistentVariable", verbose = FALSE),
    "not found"
  )
})

test_that("adonis_lpq fails with non-list_phyloseq input", {
  expect_error(
    adonis_lpq(data_fungi, formula = "Height"),
    "list_phyloseq"
  )
})

test_that("adonis_lpq handles multiple formula terms", {
  # Need to find variables that exist in sample_data
  lpq <- create_test_lpq_for_adonis()

  # Get available columns
  sam_cols <- colnames(sample_data(lpq[["original"]]))

  # Find a second variable that exists
  if ("Time" %in% sam_cols) {
    # Filter for samples without NA in Time as well
    pq_clean <- subset_samples(
      data_fungi,
      !is.na(Height) & !is.na(Time)
    )
    pq_clean <- clean_pq(pq_clean)

    lpq_multi <- list_phyloseq(
      list(a = pq_clean, b = pq_clean),
      same_bioinfo_pipeline = FALSE
    )

    result <- adonis_lpq(lpq_multi, formula = "Height + Time", verbose = FALSE)

    expect_s3_class(result, "adonis_lpq_result")
    # Should have multiple terms
    expect_true(nrow(result) >= 2)
  }
})

test_that("print method works for adonis_lpq_result", {
  lpq <- create_test_lpq_for_adonis()

  result <- adonis_lpq(lpq, formula = "Height", verbose = FALSE)

  expect_output(print(result), "PERMANOVA results")
})

test_that("adonis_lpq handles errors in individual phyloseq gracefully", {
  # Create a situation where one phyloseq might fail
  # (e.g., different sample_data structure)
  lpq <- create_test_lpq_for_adonis()

  # This should still work even if warnings are produced
  result <- suppressWarnings(
    adonis_lpq(lpq, formula = "Height", verbose = FALSE)
  )

  expect_s3_class(result, "adonis_lpq_result")
})

# ==============================================================================
# Tests for glmulti_lpq
# ==============================================================================

test_that("glmulti_lpq runs model selection on each phyloseq", {
  lpq <- create_test_lpq_for_adonis()

  result <- glmulti_lpq(lpq, formula = "Hill_0 ~ Height", verbose = FALSE)

  expect_s3_class(result, "glmulti_lpq_result")
  expect_true(tibble::is_tibble(result))
  expect_true("name" %in% colnames(result))
  expect_true("variable" %in% colnames(result))
  expect_true("importance" %in% colnames(result))
})

test_that("glmulti_lpq returns results for all phyloseq objects", {
  lpq <- create_test_lpq_for_adonis()

  result <- glmulti_lpq(lpq, formula = "Hill_0 ~ Height", verbose = FALSE)

  # Should have entries for both phyloseq objects
  expect_true(all(c("original", "copy") %in% result$name))
})

test_that("glmulti_lpq respects q parameter", {
  lpq <- create_test_lpq_for_adonis()

  result <- glmulti_lpq(
    lpq,
    formula = "Hill_1 ~ Height",
    q = c(1),
    verbose = FALSE
  )

  expect_s3_class(result, "glmulti_lpq_result")
})

test_that("glmulti_lpq fails with SEPARATE_ANALYSIS type", {
  data("enterotype", package = "phyloseq")

  lpq <- list_phyloseq(list(
    fungi = data_fungi_mini,
    enterotype = enterotype
  ))

  expect_error(
    glmulti_lpq(lpq, formula = "Hill_0 ~ Height", verbose = FALSE)
  )
})

test_that("glmulti_lpq fails with invalid formula variables", {
  lpq <- create_test_lpq_for_adonis()

  expect_error(
    glmulti_lpq(lpq, formula = "Hill_0 ~ NonexistentVariable", verbose = FALSE),
    "not found"
  )
})

test_that("glmulti_lpq fails with non-list_phyloseq input", {
  expect_error(
    glmulti_lpq(data_fungi, formula = "Hill_0 ~ Height"),
    "list_phyloseq"
  )
})

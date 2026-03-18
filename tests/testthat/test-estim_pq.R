# Tests for estim_pq.R and estim_lpq.R
# Tests for estimation statistics functions

# ==============================================================================
# Helper to create clean test data
# ==============================================================================

create_test_pq <- function() {
  pq <- subset_samples(data_fungi, !is.na(Height))
  clean_pq(pq)
}

create_test_lpq <- function() {
  pq <- create_test_pq()
  list_phyloseq(
    list(original = pq, copy = pq),
    same_bioinfo_pipeline = FALSE
  )
}

# ==============================================================================
# Tests for hill_samples_pq (internal)
# ==============================================================================

test_that("hill_samples_pq returns correct structure", {
  pq <- create_test_pq()
  res <- hill_samples_pq(pq, q = c(0, 1, 2))

  expect_true(is.data.frame(res))
  expect_equal(nrow(res), nsamples(pq))
  expect_true(all(c("Hill_0", "Hill_1", "Hill_2") %in% colnames(res)))
})

test_that("hill_samples_pq includes sample_data columns", {
  pq <- create_test_pq()
  res <- hill_samples_pq(pq)

  sam_cols <- colnames(sample_data(pq))
  expect_true(all(sam_cols %in% colnames(res)))
})

test_that("hill_samples_pq returns positive values", {
  pq <- create_test_pq()
  res <- hill_samples_pq(pq)

  expect_true(all(res$Hill_0 > 0))
  expect_true(all(res$Hill_1 > 0))
  expect_true(all(res$Hill_2 > 0))
})

test_that("hill_samples_pq respects custom q", {
  pq <- create_test_pq()
  res <- hill_samples_pq(pq, q = c(0))

  expect_true("Hill_0" %in% colnames(res))
  expect_false("Hill_1" %in% colnames(res))
  expect_false("Hill_2" %in% colnames(res))
})

# ==============================================================================
# Tests for diversity_samples_pq (internal)
# ==============================================================================

test_that("diversity_samples_pq delegates to hill when custom_fn is NULL", {
  pq <- create_test_pq()
  res <- diversity_samples_pq(pq, q = c(0, 1))

  expect_true("Hill_0" %in% colnames(res))
  expect_true("Hill_1" %in% colnames(res))
})

test_that("diversity_samples_pq works with custom function (named vector)", {
  pq <- create_test_pq()
  custom <- function(physeq) {
    sums <- sample_sums(physeq)
    stats::setNames(log(sums), sample_names(physeq))
  }

  res <- diversity_samples_pq(pq, custom_fn = custom)

  expect_true("custom_metric" %in% colnames(res))
  expect_equal(nrow(res), nsamples(pq))
})

test_that("diversity_samples_pq works with custom function (data.frame)", {
  pq <- create_test_pq()
  custom <- function(physeq) {
    data.frame(
      my_metric = log(sample_sums(physeq)),
      row.names = sample_names(physeq)
    )
  }

  res <- diversity_samples_pq(pq, custom_fn = custom)

  expect_true("my_metric" %in% colnames(res))
})

test_that("diversity_samples_pq errors on bad custom function", {
  pq <- create_test_pq()
  bad_fn <- function(physeq) "not a vector"

  expect_error(diversity_samples_pq(pq, custom_fn = bad_fn), "must return")
})

# ==============================================================================
# Tests for bootstrap_cor (internal)
# ==============================================================================

test_that("bootstrap_cor returns correct structure", {
  set.seed(42)
  x <- rnorm(50)
  y <- x + rnorm(50, sd = 0.5)
  res <- bootstrap_cor(x, y, resamples = 500)

  expect_true(is.list(res))
  expect_true(all(
    c("estimate", "ci_lower", "ci_upper", "boot_distribution") %in% names(res)
  ))
  expect_true(res$ci_lower < res$estimate)
  expect_true(res$ci_upper > res$estimate)
  expect_length(res$boot_distribution, 500)
})

# ==============================================================================
# Tests for bootstrap_lm (internal)
# ==============================================================================

test_that("bootstrap_lm returns correct structure", {
  set.seed(42)
  x <- rnorm(50)
  y <- 2 + 3 * x + rnorm(50, sd = 0.5)
  res <- bootstrap_lm(x, y, resamples = 500)

  expect_true(is.list(res))
  expect_true(all(
    c(
      "intercept",
      "slope",
      "slope_ci_lower",
      "slope_ci_upper",
      "boot_slopes"
    ) %in%
      names(res)
  ))
  expect_true(res$slope_ci_lower < res$slope)
  expect_true(res$slope_ci_upper > res$slope)
  expect_length(res$boot_slopes, 500)
})

# ==============================================================================
# Tests for estim_diff_pq
# ==============================================================================

test_that("estim_diff_pq returns correct result structure", {
  skip_if_not_installed("dabestr")
  pq <- create_test_pq()

  res <- estim_diff_pq(pq, fact = "Height", resamples = 100)

  expect_s3_class(res, "estim_diff_pq_result")
  expect_true(all(
    c("data", "dabest_objects", "plots", "summary") %in% names(res)
  ))
  expect_true(tibble::is_tibble(res$summary))
  expect_true(all(
    c(
      "metric",
      "comparison",
      "effect_size",
      "ci_lower",
      "ci_upper",
      "pvalue_permtest",
      "pvalue_welch",
      "pvalue_mann_whitney"
    ) %in%
      colnames(res$summary)
  ))
})

test_that("estim_diff_pq summary has rows for each metric", {
  skip_if_not_installed("dabestr")
  pq <- create_test_pq()

  res <- estim_diff_pq(
    pq,
    fact = "Height",
    q = c(0, 1),
    resamples = 100
  )

  metrics <- unique(res$summary$metric)
  expect_true("Hill_0" %in% metrics)
  expect_true("Hill_1" %in% metrics)
})

test_that("estim_diff_pq supports different effect types", {
  skip_if_not_installed("dabestr")
  pq <- create_test_pq()

  for (etype in c("mean_diff", "cohens_d")) {
    res <- estim_diff_pq(
      pq,
      fact = "Height",
      q = c(0),
      effect_type = etype,
      resamples = 100
    )
    expect_s3_class(res, "estim_diff_pq_result")
  }
})

test_that("estim_diff_pq handles NA removal", {
  skip_if_not_installed("dabestr")
  pq <- data_fungi # Has NAs in Height

  res <- estim_diff_pq(
    pq,
    fact = "Height",
    q = c(0),
    na_remove = TRUE,
    resamples = 100
  )
  expect_s3_class(res, "estim_diff_pq_result")
})

test_that("estim_diff_pq errors on bad input", {
  skip_if_not_installed("dabestr")
  pq <- create_test_pq()

  expect_error(estim_diff_pq(pq, fact = "nonexistent"), "not found")
  expect_error(estim_diff_pq("not_phyloseq", fact = "Height"), "phyloseq")
  expect_error(
    estim_diff_pq(pq, fact = "Height", effect_type = "bad"),
    "effect_type"
  )
})

test_that("estim_diff_pq works with custom function", {
  skip_if_not_installed("dabestr")
  pq <- create_test_pq()

  custom <- function(physeq) {
    stats::setNames(sample_sums(physeq), sample_names(physeq))
  }

  res <- estim_diff_pq(pq, fact = "Height", custom_fn = custom, resamples = 100)
  expect_s3_class(res, "estim_diff_pq_result")
  expect_true("custom_metric" %in% unique(res$summary$metric))
})

test_that("print works for estim_diff_pq_result", {
  skip_if_not_installed("dabestr")
  pq <- create_test_pq()

  res <- estim_diff_pq(pq, fact = "Height", q = c(0), resamples = 100)
  expect_output(print(res), "categorical comparison")
  expect_output(print(res), "legacy purposes")
})

# ==============================================================================
# Tests for estim_cor_pq
# ==============================================================================

test_that("estim_cor_pq returns correct result structure", {
  pq <- create_test_pq()

  # Add numeric variable
  sam <- sample_data(pq)
  sam$lib_size <- sample_sums(pq)
  sample_data(pq) <- sam

  set.seed(42)
  res <- estim_cor_pq(pq, variable = "lib_size", resamples = 100)

  expect_s3_class(res, "estim_cor_pq_result")
  expect_true(all(
    c("data", "correlations", "regressions", "plots") %in% names(res)
  ))
  expect_true(tibble::is_tibble(res$correlations))
  expect_true(tibble::is_tibble(res$regressions))
})

test_that("estim_cor_pq correlations tibble has expected columns", {
  pq <- create_test_pq()
  sam <- sample_data(pq)
  sam$lib_size <- sample_sums(pq)
  sample_data(pq) <- sam

  set.seed(42)
  res <- estim_cor_pq(
    pq,
    variable = "lib_size",
    q = c(0),
    resamples = 100
  )

  expect_true(all(
    c("metric", "estimate", "ci_lower", "ci_upper", "method", "pvalue") %in%
      colnames(res$correlations)
  ))
})

test_that("estim_cor_pq regressions tibble has expected columns", {
  pq <- create_test_pq()
  sam <- sample_data(pq)
  sam$lib_size <- sample_sums(pq)
  sample_data(pq) <- sam

  set.seed(42)
  res <- estim_cor_pq(
    pq,
    variable = "lib_size",
    q = c(0),
    resamples = 100
  )

  expect_true(all(
    c("metric", "intercept", "slope", "slope_ci_lower", "slope_ci_upper") %in%
      colnames(res$regressions)
  ))
})

test_that("estim_cor_pq supports spearman method", {
  pq <- create_test_pq()
  sam <- sample_data(pq)
  sam$lib_size <- sample_sums(pq)
  sample_data(pq) <- sam

  set.seed(42)
  res <- estim_cor_pq(
    pq,
    variable = "lib_size",
    q = c(0),
    method = "spearman",
    resamples = 100
  )

  expect_equal(res$correlations$method[1], "spearman")
})

test_that("estim_cor_pq errors on non-numeric variable", {
  pq <- create_test_pq()
  expect_error(estim_cor_pq(pq, variable = "Height"), "must be numeric")
})

test_that("estim_cor_pq errors on missing variable", {
  pq <- create_test_pq()
  expect_error(estim_cor_pq(pq, variable = "nonexistent"), "not found")
})

test_that("print works for estim_cor_pq_result", {
  pq <- create_test_pq()
  sam <- sample_data(pq)
  sam$lib_size <- sample_sums(pq)
  sample_data(pq) <- sam

  set.seed(42)
  res <- estim_cor_pq(
    pq,
    variable = "lib_size",
    q = c(0),
    resamples = 100
  )
  expect_output(print(res), "numeric correlation")
  expect_output(print(res), "legacy purposes")
})

# ==============================================================================
# Tests for estim_diff_lpq
# ==============================================================================

test_that("estim_diff_lpq returns combined summary with name column", {
  skip_if_not_installed("dabestr")
  lpq <- create_test_lpq()

  res <- estim_diff_lpq(
    lpq,
    fact = "Height",
    q = c(0),
    resamples = 100,
    verbose = FALSE
  )

  expect_s3_class(res, "estim_diff_lpq_result")
  expect_true("name" %in% colnames(res$summary))
  expect_true(all(c("original", "copy") %in% res$summary$name))
})

test_that("estim_diff_lpq rejects SEPARATE_ANALYSIS", {
  skip_if_not_installed("dabestr")
  data("enterotype", package = "phyloseq")

  lpq <- list_phyloseq(list(
    fungi = data_fungi_mini,
    enterotype = enterotype
  ))

  expect_error(
    estim_diff_lpq(lpq, fact = "Height", verbose = FALSE)
  )
})

test_that("estim_diff_lpq validates fact column", {
  skip_if_not_installed("dabestr")
  lpq <- create_test_lpq()

  expect_error(
    estim_diff_lpq(lpq, fact = "nonexistent", verbose = FALSE),
    "not found"
  )
})

test_that("estim_diff_lpq fails with non-list_phyloseq", {
  skip_if_not_installed("dabestr")
  expect_error(
    estim_diff_lpq(data_fungi, fact = "Height"),
    "list_phyloseq"
  )
})

test_that("print works for estim_diff_lpq_result", {
  skip_if_not_installed("dabestr")
  lpq <- create_test_lpq()

  res <- estim_diff_lpq(
    lpq,
    fact = "Height",
    q = c(0),
    resamples = 100,
    verbose = FALSE
  )
  expect_output(print(res), "list_phyloseq")
})

# ==============================================================================
# Tests for estim_cor_lpq
# ==============================================================================

test_that("estim_cor_lpq returns combined results with name column", {
  pq <- create_test_pq()
  sam <- sample_data(pq)
  sam$lib_size <- sample_sums(pq)
  sample_data(pq) <- sam

  lpq <- list_phyloseq(
    list(original = pq, copy = pq),
    same_bioinfo_pipeline = FALSE
  )

  set.seed(42)
  res <- estim_cor_lpq(
    lpq,
    variable = "lib_size",
    q = c(0),
    resamples = 100,
    verbose = FALSE
  )

  expect_s3_class(res, "estim_cor_lpq_result")
  expect_true("name" %in% colnames(res$correlations))
  expect_true(all(c("original", "copy") %in% res$correlations$name))
})

test_that("estim_cor_lpq validates variable column", {
  lpq <- create_test_lpq()

  expect_error(
    estim_cor_lpq(lpq, variable = "nonexistent", verbose = FALSE),
    "not found"
  )
})

test_that("estim_cor_lpq fails with non-list_phyloseq", {
  expect_error(
    estim_cor_lpq(data_fungi, variable = "lib_size"),
    "list_phyloseq"
  )
})

test_that("print works for estim_cor_lpq_result", {
  pq <- create_test_pq()
  sam <- sample_data(pq)
  sam$lib_size <- sample_sums(pq)
  sample_data(pq) <- sam

  lpq <- list_phyloseq(
    list(original = pq, copy = pq),
    same_bioinfo_pipeline = FALSE
  )

  set.seed(42)
  res <- estim_cor_lpq(
    lpq,
    variable = "lib_size",
    q = c(0),
    resamples = 100,
    verbose = FALSE
  )
  expect_output(print(res), "list_phyloseq")
})

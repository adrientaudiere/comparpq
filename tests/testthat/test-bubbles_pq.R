# Tests for bubbles_pq.R
# Tests for interactive bubble plot visualization

# ==============================================================================
# Tests for bubbles_pq
# ==============================================================================

test_that("bubbles_pq returns htmlwidget by default", {
  skip_if_not_installed("robservable")

  result <- bubbles_pq(physeq = data_fungi_mini, rank_color = "Class")

  expect_s3_class(result, "htmlwidget")
})

test_that("bubbles_pq returns dataframe when return_dataframe=TRUE", {
  result <- bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    return_dataframe = TRUE
  )

  expect_s3_class(result, "data.frame")
  expect_true("value" %in% colnames(result))
  expect_true("label" %in% colnames(result))
  expect_true("color" %in% colnames(result))
  expect_true("rank_value_color" %in% colnames(result))
  expect_true("id" %in% colnames(result))
})

test_that("bubbles_pq dataframe has correct number of rows", {
  result <- bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    return_dataframe = TRUE
  )

  expect_equal(nrow(result), ntaxa(data_fungi_mini))
})

test_that("bubbles_pq respects min_nb_seq parameter", {
  result_all <- bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    min_nb_seq = 0,
    return_dataframe = TRUE
  )

  result_filtered <- bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    min_nb_seq = 1000,
    return_dataframe = TRUE
  )

  expect_true(nrow(result_filtered) <= nrow(result_all))
})

test_that("bubbles_pq respects log1ptransform parameter", {
  result_log <- bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    log1ptransform = TRUE,
    return_dataframe = TRUE
  )

  result_nolog <- bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    log1ptransform = FALSE,
    return_dataframe = TRUE
  )

  # Log transformed values should be smaller (unless value is 0 or 1)
  # For values > 1, log1p(x) < x
  high_values <- result_nolog$value > 1
  if (any(high_values)) {
    expect_true(all(result_log$value[high_values] < result_nolog$value[high_values]))
  }
})

test_that("bubbles_pq respects randomize parameter", {
  set.seed(123)
  result1 <- bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    randomize = TRUE,
    seed = 42,
    return_dataframe = TRUE
  )

  set.seed(123)
  result2 <- bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    randomize = TRUE,
    seed = 42,
    return_dataframe = TRUE
  )

  # Same seed should give same order
  expect_equal(result1$id, result2$id)

  # Different seed should give different order
  result3 <- bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    randomize = TRUE,
    seed = 99,
    return_dataframe = TRUE
  )

  # Order should differ (most likely)
  expect_false(identical(result1$id, result3$id))
})

test_that("bubbles_pq uses taxa_names when rank_label='Taxa'", {
  result <- bubbles_pq(
    physeq = data_fungi_mini,
    rank_label = "Taxa",
    rank_color = "Class",
    return_dataframe = TRUE
  )

  expect_true(all(result$label %in% taxa_names(data_fungi_mini)))
})

test_that("bubbles_pq uses tax_table values for custom rank_label", {
  result <- bubbles_pq(
    physeq = data_fungi_mini,
    rank_label = "Genus",
    rank_color = "Class",
    return_dataframe = TRUE
  )

  expected_labels <- as.vector(data_fungi_mini@tax_table[, "Genus"])
  expect_equal(result$label, expected_labels)
})

test_that("bubbles_pq value column contains taxa_sums", {
  result <- bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    log1ptransform = FALSE,
    return_dataframe = TRUE
  )

  expected_values <- taxa_sums(data_fungi_mini)
  expect_equal(result$value, as.vector(expected_values))
})

test_that("bubbles_pq rank_value_color contains correct taxonomy", {
  result <- bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Family",
    return_dataframe = TRUE
  )

  expected_families <- as.vector(data_fungi_mini@tax_table[, "Family"])
  expect_equal(result$rank_value_color, expected_families)
})

test_that("bubbles_pq htmlwidget has correct class", {
  skip_if_not_installed("robservable")

  result <- bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class"
  )

  expect_true("robservable" %in% class(result))
})

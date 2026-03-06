test_that("gg_bubbles_pq returns a ggplot object", {
  skip_if_not_installed("packcircles")
  result <- gg_bubbles_pq(physeq = data_fungi_mini, rank_color = "Class")
  expect_s3_class(result, "ggplot")
})

test_that("gg_bubbles_pq returns dataframe when return_dataframe=TRUE", {
  result <- gg_bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    return_dataframe = TRUE
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(
    c("value", "label", "rank_value_color") %in% colnames(result)
  ))
  expect_equal(nrow(result), ntaxa(data_fungi_mini))
})

test_that("gg_bubbles_pq respects min_nb_seq", {
  result_all <- gg_bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    return_dataframe = TRUE
  )
  result_filtered <- gg_bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    min_nb_seq = 1000,
    return_dataframe = TRUE
  )
  expect_true(nrow(result_filtered) <= nrow(result_all))
})

test_that("gg_bubbles_pq respects log1ptransform", {
  result_log <- gg_bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    log1ptransform = TRUE,
    return_dataframe = TRUE
  )
  result_nolog <- gg_bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    log1ptransform = FALSE,
    return_dataframe = TRUE
  )
  high_values <- result_nolog$value > 1
  if (any(high_values)) {
    expect_true(all(
      result_log$value[high_values] < result_nolog$value[high_values]
    ))
  }
})

test_that("gg_bubbles_pq uses taxa_names when rank_label='Taxa'", {
  result <- gg_bubbles_pq(
    physeq = data_fungi_mini,
    rank_label = "Taxa",
    rank_color = "Class",
    return_dataframe = TRUE
  )
  expect_true(all(result$label %in% taxa_names(data_fungi_mini)))
})

test_that("gg_bubbles_pq uses tax_table values for custom rank_label", {
  result <- gg_bubbles_pq(
    physeq = data_fungi_mini,
    rank_label = "Genus",
    rank_color = "Class",
    return_dataframe = TRUE
  )
  expected_labels <- as.vector(data_fungi_mini@tax_table[, "Genus"])
  expect_equal(result$label, expected_labels)
})

test_that("gg_bubbles_pq facet_by adds facet column to dataframe", {
  result <- gg_bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    facet_by = "Height",
    return_dataframe = TRUE
  )
  expect_true("facet" %in% colnames(result))
  expected_levels <- unique(as.character(data_fungi_mini@sam_data[["Height"]]))
  expect_true(all(result$facet %in% expected_levels))
})

test_that("gg_bubbles_pq facet_by produces faceted ggplot", {
  skip_if_not_installed("packcircles")
  result <- gg_bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    facet_by = "Height"
  )
  expect_s3_class(result, "ggplot")
  expect_true(inherits(result$facet, "FacetWrap"))
})

test_that("gg_bubbles_pq errors on invalid facet_by", {
  expect_error(
    gg_bubbles_pq(
      physeq = data_fungi_mini,
      rank_color = "Class",
      facet_by = "nonexistent_column"
    ),
    "not found in sample_data"
  )
})

test_that("gg_bubbles_pq square layout returns ggplot", {
  skip_if_not_installed("packcircles")
  result <- gg_bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    layout = "square"
  )
  expect_s3_class(result, "ggplot")
})

test_that("gg_bubbles_pq square layout with facet returns ggplot", {
  skip_if_not_installed("packcircles")
  result <- gg_bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    layout = "square",
    facet_by = "Height"
  )
  expect_s3_class(result, "ggplot")
  expect_true(inherits(result$facet, "FacetWrap"))
})

test_that("gg_bubbles_pq rank_contour colors borders by tax_table column", {
  skip_if_not_installed("packcircles")
  result <- gg_bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    rank_contour = "Order"
  )
  expect_s3_class(result, "ggplot")
  layer_data <- ggplot2::layer_data(result, 1)
  expect_true("colour" %in% colnames(layer_data))
  expect_true(length(unique(layer_data$colour)) > 1)
})

test_that("gg_bubbles_pq rank_contour adds column to return_dataframe", {
  result <- gg_bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    rank_contour = "Order",
    return_dataframe = TRUE
  )
  expect_true("rank_value_contour" %in% colnames(result))
  expected <- as.vector(data_fungi_mini@tax_table[, "Order"])
  expect_equal(result$rank_value_contour, expected)
})

test_that("gg_bubbles_pq with show_labels=FALSE still returns ggplot", {
  skip_if_not_installed("packcircles")
  result <- gg_bubbles_pq(
    physeq = data_fungi_mini,
    rank_color = "Class",
    show_labels = FALSE
  )
  expect_s3_class(result, "ggplot")
})

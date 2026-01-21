# Tests for formattable_lpq.R
# Tests for visualization functions for list_phyloseq objects

# ==============================================================================
# Tests for formattable_lpq
# ==============================================================================

test_that("formattable_lpq creates formattable object", {
  skip_if_not_installed("formattable")

  lpq <- list_phyloseq(list(fungi = data_fungi, mini = data_fungi_mini))

  result <- formattable_lpq(lpq)

  expect_s3_class(result, "formattable")
})

test_that("formattable_lpq respects columns parameter", {
  skip_if_not_installed("formattable")

  lpq <- list_phyloseq(list(fungi = data_fungi, mini = data_fungi_mini))

  result <- formattable_lpq(lpq, columns = c("name", "n_samples", "n_taxa"))

  # Check that result has correct columns
  df <- as.data.frame(result)
  expect_equal(ncol(df), 3)
  expect_true(all(c("name", "n_samples", "n_taxa") %in% colnames(df)))
})

test_that("formattable_lpq respects void_style parameter", {
  skip_if_not_installed("formattable")

  lpq <- list_phyloseq(list(fungi = data_fungi))

  result <- formattable_lpq(lpq, void_style = TRUE)

  expect_s3_class(result, "formattable")
})

test_that("formattable_lpq respects log10_transform parameter", {
  skip_if_not_installed("formattable")

  lpq <- list_phyloseq(list(fungi = data_fungi))

  result_log <- formattable_lpq(lpq, log10_transform = TRUE)
  result_nolog <- formattable_lpq(lpq, log10_transform = FALSE)

  # Both should be formattable objects
  expect_s3_class(result_log, "formattable")
  expect_s3_class(result_nolog, "formattable")
})

test_that("formattable_lpq accepts custom bar_colors", {
  skip_if_not_installed("formattable")

  lpq <- list_phyloseq(list(fungi = data_fungi))

  result <- formattable_lpq(lpq, bar_colors = list(
    n_samples = "red",
    n_taxa = "blue"
  ))

  expect_s3_class(result, "formattable")
})

test_that("formattable_lpq fails with non-list_phyloseq input", {
  skip_if_not_installed("formattable")

  expect_error(
    formattable_lpq(data_fungi),
    "list_phyloseq"
  )
})

# ==============================================================================
# Tests for formattable_lpq_full
# ==============================================================================

test_that("formattable_lpq_full returns list with summary and comparison", {
  skip_if_not_installed("formattable")

  lpq <- list_phyloseq(list(fungi = data_fungi, mini = data_fungi_mini))

  result <- formattable_lpq_full(lpq)

  expect_type(result, "list")
  expect_true("summary" %in% names(result))
  expect_true("comparison" %in% names(result))
})

test_that("formattable_lpq_full respects show_summary parameter", {
  skip_if_not_installed("formattable")

  lpq <- list_phyloseq(list(fungi = data_fungi))

  result <- formattable_lpq_full(lpq, show_summary = FALSE, show_comparison = TRUE)

  expect_s3_class(result, "formattable")
})

test_that("formattable_lpq_full respects show_comparison parameter", {
  skip_if_not_installed("formattable")

  lpq <- list_phyloseq(list(fungi = data_fungi))

  result <- formattable_lpq_full(lpq, show_summary = TRUE, show_comparison = FALSE)

  expect_s3_class(result, "formattable")
})

# ==============================================================================
# Tests for factor_formatter
# ==============================================================================

test_that("factor_formatter handles factor input", {
  skip_if_not_installed("formattable")

  x <- factor(c("A", "B", "A", "C"))

  result <- factor_formatter(x)

  expect_type(result, "closure")
})
test_that("factor_formatter handles character input", {
  skip_if_not_installed("formattable")

  x <- c("A", "B", "A", "C")

  result <- factor_formatter(x)

  expect_type(result, "closure")
})

test_that("factor_formatter returns input for non-factor/character", {
  skip_if_not_installed("formattable")

  x <- c(1, 2, 3)

  result <- factor_formatter(x)

  expect_equal(result, x)
})

# ==============================================================================
# Tests for upset_lpq
# ==============================================================================

test_that("upset_lpq creates venn diagram for small number of sets", {
  skip_if_not_installed("ggVennDiagram")

  lpq <- list_phyloseq(list(fungi = data_fungi, mini = data_fungi_mini))

  # Should auto-select venn for 2 sets
  result <- upset_lpq(lpq, tax_rank = "Family")

  expect_s3_class(result, "ggplot")
})

test_that("upset_lpq creates upset plot when plot_type='upset'", {
  skip_if_not_installed("ComplexUpset")

  lpq <- list_phyloseq(list(fungi = data_fungi, mini = data_fungi_mini))

  result <- upset_lpq(lpq, tax_rank = "Family", plot_type = "upset")

  # ComplexUpset returns a patchwork object
  expect_true(inherits(result, "patchwork") || inherits(result, "ggplot"))
})

test_that("upset_lpq creates venn diagram when plot_type='venn'", {
  skip_if_not_installed("ggVennDiagram")

  lpq <- list_phyloseq(list(fungi = data_fungi, mini = data_fungi_mini))

  result <- upset_lpq(lpq, tax_rank = "Family", plot_type = "venn")

  expect_s3_class(result, "ggplot")
})

test_that("upset_lpq respects remove_na parameter", {
  skip_if_not_installed("ggVennDiagram")

  lpq <- list_phyloseq(list(fungi = data_fungi, mini = data_fungi_mini))

  result_with_na <- upset_lpq(lpq, tax_rank = "Family", remove_na = FALSE)
  result_without_na <- upset_lpq(lpq, tax_rank = "Family", remove_na = TRUE)

  # Both should be ggplot objects
  expect_s3_class(result_with_na, "ggplot")
  expect_s3_class(result_without_na, "ggplot")
})

test_that("upset_lpq fails with invalid tax_rank", {
  lpq <- list_phyloseq(list(fungi = data_fungi))

  expect_error(
    upset_lpq(lpq, tax_rank = "NonexistentRank"),
    "not found in tax_table"
  )
})

test_that("upset_lpq fails with non-list_phyloseq input", {
  expect_error(
    upset_lpq(data_fungi, tax_rank = "Family"),
    "list_phyloseq"
  )
})

test_that("upset_lpq handles multiple phyloseq objects", {
  skip_if_not_installed("ComplexUpset")

  set.seed(123)
  rarefied <- rarefy_even_depth(data_fungi, sample.size = 1000, verbose = FALSE)

  lpq <- list_phyloseq(list(
    fungi = data_fungi,
    mini = data_fungi_mini,
    rarefied = rarefied
  ))

  # Should auto-select venn for 3 sets
  result <- upset_lpq(lpq, tax_rank = "Family")

  expect_true(inherits(result, "ggplot") || inherits(result, "patchwork"))
})

test_that("upset_lpq warns for venn with >7 sets", {
  skip_if_not_installed("ggVennDiagram")

  # Create list with many phyloseq objects
  pq_list <- list()
  for (i in 1:8) {
    pq_list[[paste0("pq_", i)]] <- data_fungi_mini
  }
  lpq <- list_phyloseq(pq_list)

  expect_warning(
    upset_lpq(lpq, tax_rank = "Family", plot_type = "venn"),
    "more than 7 sets"
  )
})

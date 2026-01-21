# Tests for compare_taxo_plot.R
# Tests for visualization functions for taxonomic comparisons

# ==============================================================================
# Tests for tc_bar
# ==============================================================================

test_that("tc_bar creates ggplot object", {
  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000)

  result <- tc_bar(pq, rank_1 = 5, rank_2 = 13, color_rank = 3)

  expect_s3_class(result, "ggplot")
})

test_that("tc_bar works with character rank names", {
  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000)

  result <- tc_bar(
    pq,
    rank_1 = "Genus",
    rank_2 = "Genus__eukaryome_Glomero",
    color_rank = "Family"
  )

  expect_s3_class(result, "ggplot")
})

test_that("tc_bar respects log10trans parameter", {
  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000)

  result_log <- tc_bar(pq, rank_1 = 5, rank_2 = 13, color_rank = 3, log10trans = TRUE)
  result_nolog <- tc_bar(pq, rank_1 = 5, rank_2 = 13, color_rank = 3, log10trans = FALSE)

  expect_s3_class(result_log, "ggplot")
  expect_s3_class(result_nolog, "ggplot")
})

test_that("tc_bar respects point_size and point_alpha", {
  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000)

  result <- tc_bar(
    pq,
    rank_1 = 5,
    rank_2 = 13,
    color_rank = 3,
    point_size = 0.5,
    point_alpha = 0.5
  )

  expect_s3_class(result, "ggplot")
})

# ==============================================================================
# Tests for tc_points_matrix
# ==============================================================================

test_that("tc_points_matrix creates ggplot object", {
  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000)

  result <- tc_points_matrix(pq, "Order", "Order__eukaryome_Glomero")

  expect_s3_class(result, "ggplot")
})

test_that("tc_points_matrix works with integer ranks", {
  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000)

  result <- tc_points_matrix(pq, 6, 14)

  expect_s3_class(result, "ggplot")
})

test_that("tc_points_matrix respects stat_across_sample parameter", {
  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000)

  result_sum <- tc_points_matrix(pq, 4, 12, stat_across_sample = "sum")
  result_mean <- tc_points_matrix(pq, 4, 12, stat_across_sample = "mean")

  expect_s3_class(result_sum, "ggplot")
  expect_s3_class(result_mean, "ggplot")
})

test_that("tc_points_matrix fails with invalid stat_across_sample", {
  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000)

  expect_error(
    tc_points_matrix(pq, 4, 12, stat_across_sample = "invalid"),
    "must be set to mean or sum"
  )
})

test_that("tc_points_matrix respects custom colors", {
  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000)

  result <- tc_points_matrix(
    pq,
    "Order",
    "Order__eukaryome_Glomero",
    color_1 = "#ff0000",
    color_2 = "#0000ff"
  )

  expect_s3_class(result, "ggplot")
})

# ==============================================================================
# Tests for rainplot_taxo_na
# ==============================================================================

test_that("rainplot_taxo_na creates ggplot object", {
  skip_if_not_installed("ggrain")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000)

  result <- rainplot_taxo_na(
    pq,
    ranks = c("Family", "Family__eukaryome_Glomero")
  )

  expect_s3_class(result, "ggplot")
})

test_that("rainplot_taxo_na works with integer ranks", {
  skip_if_not_installed("ggrain")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000)

  result <- rainplot_taxo_na(pq, ranks = c(4, 12))

  expect_s3_class(result, "ggplot")
})

test_that("rainplot_taxo_na uses all ranks when ranks is NULL", {
  skip_if_not_installed("ggrain")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 10000)

  result <- rainplot_taxo_na(pq, ranks = NULL)

  expect_s3_class(result, "ggplot")
})

test_that("rainplot_taxo_na respects min_nb_seq parameter", {
  skip_if_not_installed("ggrain")

  pq <- Glom_otu

  result <- rainplot_taxo_na(pq, ranks = c(4, 12), min_nb_seq = 100)

  expect_s3_class(result, "ggplot")
})

test_that("rainplot_taxo_na respects sample_colored parameter", {
  skip_if_not_installed("ggrain")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000)

  result <- rainplot_taxo_na(
    pq,
    ranks = c(4, 12),
    sample_colored = TRUE
  )

  expect_s3_class(result, "ggplot")
})

test_that("rainplot_taxo_na respects sample_linked parameter", {
  skip_if_not_installed("ggrain")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000)

  result <- rainplot_taxo_na(
    pq,
    ranks = c(4, 12),
    sample_linked = TRUE
  )

  expect_s3_class(result, "ggplot")
})

# ==============================================================================
# Tests for tc_circle
# ==============================================================================

test_that("tc_circle runs without error", {
  skip_if_not_installed("circlize")

  expect_no_error(
    tc_circle(
      Glom_otu,
      "Genus__eukaryome_Glomero",
      "Genus",
      suffix_1 = "_Euk",
      suffix_2 = "_Marjaam"
    )
  )
})

test_that("tc_circle works with different ranks", {
  skip_if_not_installed("circlize")

  expect_no_error(
    tc_circle(
      Glom_otu,
      "Family__eukaryome_Glomero",
      "Family",
      suffix_1 = "_Euk",
      suffix_2 = "_Marjaam"
    )
  )
})

test_that("tc_circle respects custom suffixes", {
  skip_if_not_installed("circlize")

  expect_no_error(
    tc_circle(
      Glom_otu,
      "Genus__eukaryome_Glomero",
      "Genus",
      suffix_1 = "_A",
      suffix_2 = "_B"
    )
  )
})

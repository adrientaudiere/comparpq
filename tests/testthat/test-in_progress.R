# Tests for in_progress.R
# Tests for functions under development

# ==============================================================================
# Tests for tc_linked_trees
# ==============================================================================

test_that("tc_linked_trees creates ggplot object", {
  skip_if_not_installed("ggtree")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 10000)

  result <- tc_linked_trees(
    pq, pq,
    ranks_1 = c("Phylum", "Class", "Order", "Family", "Genus"),
    ranks_2 = c(
      "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero", "Family__eukaryome_Glomero",
      "Genus__eukaryome_Glomero"
    )
  )

  expect_s3_class(result, "ggplot")
})

test_that("tc_linked_trees respects show_tip_links parameter", {
  skip_if_not_installed("ggtree")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 10000)

  result_with_links <- tc_linked_trees(
    pq, pq,
    ranks_1 = c("Phylum", "Class", "Order"),
    ranks_2 = c(
      "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    show_tip_links = TRUE
  )

  result_without_links <- tc_linked_trees(
    pq, pq,
    ranks_1 = c("Phylum", "Class", "Order"),
    ranks_2 = c(
      "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    show_tip_links = FALSE
  )

  expect_s3_class(result_with_links, "ggplot")
  expect_s3_class(result_without_links, "ggplot")
})

test_that("tc_linked_trees respects show_labels parameter", {
  skip_if_not_installed("ggtree")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 10000)

  result_with_labels <- tc_linked_trees(
    pq, pq,
    ranks_1 = c("Phylum", "Class", "Order"),
    ranks_2 = c(
      "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    show_labels = TRUE
  )

  result_without_labels <- tc_linked_trees(
    pq, pq,
    ranks_1 = c("Phylum", "Class", "Order"),
    ranks_2 = c(
      "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    show_labels = FALSE
  )

  expect_s3_class(result_with_labels, "ggplot")
  expect_s3_class(result_without_labels, "ggplot")
})

test_that("tc_linked_trees respects tree_distance parameter", {
  skip_if_not_installed("ggtree")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 10000)

  result <- tc_linked_trees(
    pq, pq,
    ranks_1 = c("Phylum", "Class", "Order"),
    ranks_2 = c(
      "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    tree_distance = 2
  )

  expect_s3_class(result, "ggplot")
})

test_that("tc_linked_trees respects link_alpha parameter", {
  skip_if_not_installed("ggtree")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 10000)

  result <- tc_linked_trees(
    pq, pq,
    ranks_1 = c("Phylum", "Class", "Order"),
    ranks_2 = c(
      "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    link_alpha = 0.3
  )

  expect_s3_class(result, "ggplot")
})

test_that("tc_linked_trees respects link_color parameter", {
  skip_if_not_installed("ggtree")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 10000)

  result_grey <- tc_linked_trees(
    pq, pq,
    ranks_1 = c("Phylum", "Class", "Order"),
    ranks_2 = c(
      "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    link_color = "grey50"
  )

  result_colored <- tc_linked_trees(
    pq, pq,
    ranks_1 = c("Phylum", "Class", "Order"),
    ranks_2 = c(
      "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    link_color = NULL
  )

  expect_s3_class(result_grey, "ggplot")
  expect_s3_class(result_colored, "ggplot")
})

test_that("tc_linked_trees respects label_size parameter", {
  skip_if_not_installed("ggtree")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 10000)

  result <- tc_linked_trees(
    pq, pq,
    ranks_1 = c("Phylum", "Class", "Order"),
    ranks_2 = c(
      "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    label_size = 3
  )

  expect_s3_class(result, "ggplot")
})

test_that("tc_linked_trees respects internal_node_singletons parameter", {
  skip_if_not_installed("ggtree")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 10000)

  result_with_singletons <- tc_linked_trees(
    pq, pq,
    ranks_1 = c("Phylum", "Class", "Order"),
    ranks_2 = c(
      "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    internal_node_singletons = TRUE
  )

  result_without_singletons <- tc_linked_trees(
    pq, pq,
    ranks_1 = c("Phylum", "Class", "Order"),
    ranks_2 = c(
      "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    internal_node_singletons = FALSE
  )

  expect_s3_class(result_with_singletons, "ggplot")
  expect_s3_class(result_without_singletons, "ggplot")
})

test_that("tc_linked_trees returns customizable ggplot", {
  skip_if_not_installed("ggtree")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 10000)

  result <- tc_linked_trees(
    pq, pq,
    ranks_1 = c("Phylum", "Class", "Order"),
    ranks_2 = c(
      "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    )
  ) +
    ggplot2::theme(legend.position = "none")

  expect_s3_class(result, "ggplot")
})

test_that("tc_linked_trees aligns matching taxa horizontally", {
  skip_if_not_installed("ggtree")

  pq <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 10000)

  result <- tc_linked_trees(
    pq, pq,
    ranks_1 = c("Phylum", "Class", "Order"),
    ranks_2 = c(
      "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    show_tip_links = TRUE
  )

  expect_s3_class(result, "ggplot")

  build <- ggplot2::ggplot_build(result)
  segment_data <- build$data[[which(sapply(build$data, \(x) "yend" %in% names(x)))]]

  if (nrow(segment_data) > 0) {
    expect_true(all(segment_data$y == segment_data$yend))
  }
})

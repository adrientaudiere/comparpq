# Tests for tc_linked_trees.R
# Tests for taxonomic tree comparison functions

# ==============================================================================
# Tests for tc_congruence_metrics
# ==============================================================================

test_that("tc_congruence_metrics returns correct structure", {
  result <- tc_congruence_metrics(
    Glom_otu,
    ranks_1 = c("Kingdom", "Class", "Order", "Family"),
    ranks_2 = c(
      "Kingdom__eukaryome_Glomero",
      "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero",
      "Family__eukaryome_Glomero"
    )
  )

  expect_type(result, "list")
  expect_true("summary" %in% names(result))
  expect_true("total_congruent" %in% names(result))
  expect_true("partial_1_deeper" %in% names(result))
  expect_true("partial_2_deeper" %in% names(result))
  expect_true("incongruent_leaves" %in% names(result))
  expect_true("incongruent_nodes" %in% names(result))
  expect_true("details" %in% names(result))

  expect_s3_class(result$summary, "data.frame")
  expect_s3_class(result$details, "data.frame")
})

test_that("tc_congruence_metrics summary has correct columns", {
  result <- tc_congruence_metrics(
    Glom_otu,
    ranks_1 = c("Kingdom", "Class", "Order", "Family"),
    ranks_2 = c(
      "Kingdom__eukaryome_Glomero",
      "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero",
      "Family__eukaryome_Glomero"
    )
  )

  expect_true("category" %in% names(result$summary))
  expect_true("count" %in% names(result$summary))
  expect_true("percentage" %in% names(result$summary))
})

test_that("tc_congruence_metrics counts are non-negative", {
  result <- tc_congruence_metrics(
    Glom_otu,
    ranks_1 = c("Kingdom", "Class", "Order", "Family"),
    ranks_2 = c(
      "Kingdom__eukaryome_Glomero",
      "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero",
      "Family__eukaryome_Glomero"
    )
  )

  expect_true(all(result$summary$count >= 0))
  expect_true(all(result$summary$percentage >= 0))
  expect_true(all(result$summary$percentage <= 100))
})

test_that("tc_congruence_metrics total counts match number of taxa", {
  result <- tc_congruence_metrics(
    Glom_otu,
    ranks_1 = c("Kingdom", "Class", "Order", "Family"),
    ranks_2 = c(
      "Kingdom__eukaryome_Glomero",
      "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero",
      "Family__eukaryome_Glomero"
    )
  )

  n_taxa <- phyloseq::ntaxa(Glom_otu)
  total_in_details <- nrow(result$details)

  # All common taxa should be accounted for in details

  expect_equal(total_in_details, n_taxa)
})

test_that("tc_congruence_metrics details has correct columns", {
  result <- tc_congruence_metrics(
    Glom_otu,
    ranks_1 = c("Kingdom", "Class", "Order", "Family"),
    ranks_2 = c(
      "Kingdom__eukaryome_Glomero",
      "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero",
      "Family__eukaryome_Glomero"
    )
  )

  expect_true("taxon" %in% names(result$details))
  expect_true("depth_1" %in% names(result$details))
  expect_true("depth_2" %in% names(result$details))
  expect_true("category" %in% names(result$details))
})

test_that("tc_congruence_metrics works with NULL physeq_2", {
  result <- tc_congruence_metrics(
    Glom_otu,
    physeq_2 = NULL,
    ranks_1 = c("Kingdom", "Class", "Order", "Family"),
    ranks_2 = c(
      "Kingdom__eukaryome_Glomero",
      "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero",
      "Family__eukaryome_Glomero"
    )
  )

  expect_type(result, "list")
  expect_equal(length(result$only_in_1), 0)
  expect_equal(length(result$only_in_2), 0)
})

test_that("tc_congruence_metrics categories are mutually exclusive", {
  result <- tc_congruence_metrics(
    Glom_otu,
    ranks_1 = c("Kingdom", "Class", "Order", "Family"),
    ranks_2 = c(
      "Kingdom__eukaryome_Glomero",
      "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero",
      "Family__eukaryome_Glomero"
    )
  )

  all_categorized <- c(
    result$classified_only_1,
    result$classified_only_2,
    result$unclassified_both,
    result$congruent,
    result$partial_1_deeper,
    result$partial_2_deeper,
    result$incongruent_leaves,
    result$incongruent_nodes
  )

  # No duplicates in categorization
  expect_equal(length(all_categorized), length(unique(all_categorized)))
})

test_that("tc_congruence_metrics warns with no common taxa", {
  # Create a subset with different taxa names
  physeq_subset <- phyloseq::prune_taxa(taxa_names(Glom_otu)[1:5], Glom_otu)
  taxa_names(physeq_subset) <- paste0("new_", taxa_names(physeq_subset))

  expect_warning(
    tc_congruence_metrics(
      Glom_otu,
      physeq_2 = physeq_subset,
      ranks_1 = c("Kingdom", "Class"),
      ranks_2 = c("Kingdom", "Class")
    ),
    "No common taxa names"
  )
})

test_that("tc_congruence_metrics works with some common taxa", {
  # Create a subset with different taxa names
  physeq_subset <- phyloseq::prune_taxa(taxa_names(Glom_otu)[1:5], Glom_otu)
  taxa_names(physeq_subset)[1:2] <- paste0(
    "new_",
    taxa_names(physeq_subset)[1:2]
  )

  result <- tc_congruence_metrics(
    Glom_otu,
    physeq_2 = physeq_subset,
    ranks_1 = c("Kingdom", "Class"),
    ranks_2 = c("Kingdom", "Class")
  )
  expect_type(result, "list")

  expect_equal(
    result$summary |> filter(category == "total_congruent") |> pull(count),
    3
  )
})


# ==============================================================================
# Tests for tc_linked_trees
# ==============================================================================

test_that("tc_linked_trees returns a ggplot object", {
  skip_if_not_installed("ggtree")

  p <- tc_linked_trees(
    Glom_otu,
    ranks_1 = c("Kingdom", "Class", "Order"),
    ranks_2 = c(
      "Kingdom__eukaryome_Glomero",
      "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    )
  )

  expect_s3_class(p, "ggplot")
})

test_that("tc_linked_trees works with link_by_taxa = TRUE", {
  skip_if_not_installed("ggtree")

  p <- tc_linked_trees(
    Glom_otu,
    ranks_1 = c("Kingdom", "Class", "Order"),
    ranks_2 = c(
      "Kingdom__eukaryome_Glomero",
      "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    link_by_taxa = TRUE
  )

  expect_s3_class(p, "ggplot")
})

test_that("tc_linked_trees works with show_tip_links = FALSE", {
  skip_if_not_installed("ggtree")

  p <- tc_linked_trees(
    Glom_otu,
    ranks_1 = c("Kingdom", "Class", "Order"),
    ranks_2 = c(
      "Kingdom__eukaryome_Glomero",
      "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    show_tip_links = FALSE
  )

  expect_s3_class(p, "ggplot")
})

test_that("tc_linked_trees works with show_labels = FALSE", {
  skip_if_not_installed("ggtree")

  p <- tc_linked_trees(
    Glom_otu,
    ranks_1 = c("Kingdom", "Class", "Order"),
    ranks_2 = c(
      "Kingdom__eukaryome_Glomero",
      "Class__eukaryome_Glomero",
      "Order__eukaryome_Glomero"
    ),
    show_labels = FALSE
  )

  expect_s3_class(p, "ggplot")
})

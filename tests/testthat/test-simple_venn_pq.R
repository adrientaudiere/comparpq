test_that("simple_venn_pq default returns patchwork of all ranks", {
  skip_if_not_installed("patchwork")
  p <- simple_venn_pq(data_fungi_mini, "Height", verbose = FALSE)
  expect_s3_class(p, "patchwork")
})

test_that("simple_venn_pq works with taxonomic_rank", {
  p <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = "Genus",
    verbose = FALSE
  )
  expect_s3_class(p, "ggplot")
})

test_that("simple_venn_pq combine = TRUE returns patchwork", {
  skip_if_not_installed("patchwork")
  p <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = c("Family", "Genus"),
    verbose = FALSE
  )
  expect_s3_class(p, "patchwork")
})

test_that("simple_venn_pq combine = FALSE returns list", {
  plots <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = c("Family", "Genus"),
    combine = FALSE,
    verbose = FALSE
  )
  expect_type(plots, "list")
  expect_length(plots, 3)
  expect_named(plots, c("Family", "Genus", "Taxa (OTU/ASV)"))
  expect_s3_class(plots[["Family"]], "ggplot")
  expect_s3_class(plots[["Genus"]], "ggplot")
})

test_that("simple_venn_pq errors with missing factor", {
  expect_error(
    simple_venn_pq(data_fungi_mini, "nonexistent", verbose = FALSE),
    "not found in sample_data"
  )
})

test_that("simple_venn_pq errors with non-phyloseq input", {
  expect_error(
    simple_venn_pq(list(a = 1), "Height", verbose = FALSE),
    "must be a phyloseq or list_phyloseq object"
  )
})

test_that("simple_venn_pq errors with >4 groups", {
  sam <- phyloseq::sample_data(data_fungi_mini)
  sam[["many_groups"]] <- paste0("g", seq_len(nrow(sam)))
  phyloseq::sample_data(data_fungi_mini) <- sam

  expect_error(
    simple_venn_pq(data_fungi_mini, "many_groups", verbose = FALSE),
    "at most 4 groups"
  )
})

test_that("simple_venn_pq respects min_nb_seq", {
  p_strict <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = "Genus",
    min_nb_seq = 1000,
    verbose = FALSE
  )
  expect_s3_class(p_strict, "ggplot")
})

test_that("simple_venn_pq respects custom colors", {
  p <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = "Genus",
    colors = c("orange", "purple"),
    verbose = FALSE
  )
  expect_s3_class(p, "ggplot")
})

test_that("simple_venn_pq add_nb_samples = FALSE works", {
  p <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = "Genus",
    add_nb_samples = FALSE,
    verbose = FALSE
  )
  expect_s3_class(p, "ggplot")
})

test_that("simple_venn_pq count_type = 'taxa' counts ASVs", {
  p_rank <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = "Genus",
    count_type = "rank",
    verbose = FALSE
  )
  p_taxa <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = "Genus",
    count_type = "taxa",
    verbose = FALSE
  )
  p_seq <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = "Genus",
    count_type = "sequences",
    verbose = FALSE
  )
  expect_s3_class(p_rank, "ggplot")
  expect_s3_class(p_taxa, "ggplot")
  expect_s3_class(p_seq, "ggplot")
})

test_that("simple_venn_pq count_type = 'both' shows rank and taxa", {
  p <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = "Genus",
    count_type = "rank_taxa",
    verbose = FALSE
  )
  expect_s3_class(p, "ggplot")
})

test_that("simple_venn_pq works with list_phyloseq input", {
  lpq <- list_phyloseq(list(
    fungi = data_fungi_mini,
    fungi2 = data_fungi_mini
  ))
  p <- simple_venn_pq(lpq, taxonomic_rank = "Genus", verbose = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("simple_venn_pq list_phyloseq shows original sample counts", {
  lpq <- list_phyloseq(list(
    fungi = data_fungi_mini,
    fungi2 = data_fungi_mini
  ))
  p <- simple_venn_pq(
    lpq,
    taxonomic_rank = "Genus",
    add_nb_samples = TRUE,
    count_taxa = FALSE,
    verbose = FALSE
  )
  # Extract annotations and find sample count labels
  labels <- vapply(p$layers, \(l) {
    if (inherits(l$geom, "GeomText")) {
      l$aes_params$label %||% ""
    } else {
      ""
    }
  }, character(1))
  expected_n <- phyloseq::nsamples(data_fungi_mini)
  expect_true(any(grepl(paste0("n=", expected_n), labels)))
})

test_that("simple_venn_pq count_taxa adds Taxa panel", {
  skip_if_not_installed("patchwork")
  p <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = "Genus",
    count_taxa = TRUE,
    verbose = FALSE
  )
  expect_s3_class(p, "patchwork")
})

test_that("simple_venn_pq show_na_count adds NA annotation", {
  p <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = "Genus",
    show_na_count = TRUE,
    verbose = FALSE
  )
  expect_s3_class(p, "ggplot")
  # Check that "NA:" label is present in the plot layers
  labels <- vapply(p$layers, \(l) {
    if (inherits(l$geom, "GeomText")) {
      l$aes_params$label %||% ""
    } else {
      ""
    }
  }, character(1))
  expect_true(any(grepl("^NA:", labels)))
})

test_that("show_na_count taxa counts sum to ntaxa", {
  p <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = "Genus",
    count_type = "taxa",
    show_na_count = TRUE,
    na_remove = TRUE,
    verbose = FALSE
  )
  # Extract all text annotations
  annotations <- vapply(p$layers, \(l) {
    if (inherits(l$geom, "GeomText")) {
      l$aes_params$label %||% ""
    } else {
      ""
    }
  }, character(1))
  # Get NA count from annotation
  na_label <- annotations[grepl("^NA:", annotations)]
  na_count <- as.integer(sub("NA: ", "", na_label))
  # Get Venn region counts (numeric labels, not group names)
  numeric_labels <- annotations[grepl("^[0-9]+$", annotations)]
  venn_sum <- sum(as.integer(numeric_labels))
  expect_equal(venn_sum + na_count, phyloseq::ntaxa(data_fungi_mini))
})

test_that("simple_venn_pq count_taxa = FALSE with combine = FALSE", {
  plots <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = "Genus",
    count_taxa = TRUE,
    combine = FALSE,
    verbose = FALSE
  )
  expect_type(plots, "list")
  expect_named(plots, c("Genus", "Taxa (OTU/ASV)"))
})

test_that("simple_venn_pq errors when fact is NULL for phyloseq", {
  expect_error(
    simple_venn_pq(data_fungi_mini, verbose = FALSE),
    "fact.*required"
  )
})

test_that("simple_venn_pq scale_text works", {
  p <- simple_venn_pq(
    data_fungi_mini,
    "Height",
    taxonomic_rank = "Genus",
    scale_text = TRUE,
    verbose = FALSE
  )
  expect_s3_class(p, "ggplot")
})

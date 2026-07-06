skip_on_cran()
library(MiscMetabar)
data(data_fungi_mini)

small_pq <- function() {
  pq <- prune_samples(sample_names(data_fungi_mini)[1:4], data_fungi_mini)
  prune_taxa(taxa_names(pq)[1:20], pq)
}

test_that("merge_clust_lpq keeps all samples and adds a source column", {
  skip_if_not(MiscMetabar::is_vsearch_installed())
  pq <- small_pq()
  lpq <- list_phyloseq(list(run_a = pq, run_b = pq))

  merged <- merge_clust_lpq(lpq, id = 0.97, verbose = FALSE)

  expect_s4_class(merged, "phyloseq")
  expect_equal(phyloseq::nsamples(merged), 2 * phyloseq::nsamples(pq))
  expect_true("source_name" %in% colnames(phyloseq::sample_data(merged)))
  expect_setequal(
    unique(as.character(phyloseq::sample_data(merged)$source_name)),
    c("run_a", "run_b")
  )
})

test_that("merge_clust_lpq suffixes colliding sample names by parent object", {
  skip_if_not(MiscMetabar::is_vsearch_installed())
  pq <- small_pq()
  merged <- merge_clust_lpq(
    list(run_a = pq, run_b = pq),
    id = 0.97,
    verbose = FALSE
  )
  smp <- phyloseq::sample_names(merged)
  expect_equal(length(smp), length(unique(smp)))
  expect_true(any(grepl("_run_a$", smp)))
  expect_true(any(grepl("_run_b$", smp)))
})

test_that("merge_clust_lpq clusters duplicated sequences together", {
  skip_if_not(MiscMetabar::is_vsearch_installed())
  pq <- small_pq()
  merged <- merge_clust_lpq(
    list(run_a = pq, run_b = pq),
    id = 0.97,
    verbose = FALSE
  )
  # Two identical objects -> merged taxa count should not exceed the pooled
  # count and should be at most the single-object taxa count.
  expect_lte(phyloseq::ntaxa(merged), phyloseq::ntaxa(pq))
})

test_that("merge_clust_lpq errors on too few objects or missing refseq", {
  pq <- small_pq()
  expect_error(
    merge_clust_lpq(list(run_a = pq), verbose = FALSE),
    "at least two"
  )

  pq_no_seq <- pq
  pq_no_seq@refseq <- NULL
  expect_error(
    merge_clust_lpq(list(a = pq_no_seq, b = pq_no_seq), verbose = FALSE),
    "refseq"
  )
})

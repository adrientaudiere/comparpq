test_that("compare_refseq_pq works with identical objects", {
  res <- compare_refseq_pq(data_fungi_mini, data_fungi_mini, verbose = FALSE)

  expect_s3_class(res, "compare_refseq")
  expect_equal(length(res$shared_names), ntaxa(data_fungi_mini))
  expect_equal(length(res$unique_names_1), 0)
  expect_equal(length(res$unique_names_2), 0)
  expect_equal(length(res$unique_seqs_1), 0)
  expect_equal(length(res$unique_seqs_2), 0)
  expect_equal(nrow(res$same_name_diff_seq), 0)
  expect_equal(nrow(res$same_seq_diff_name), 0)
  expect_true(is.na(res$mean_nn_kmer_dist_1))
  expect_true(is.na(res$mean_nn_kmer_dist_2))
})

test_that("compare_refseq_pq detects unique taxa and computes distances", {
  taxa_keep <- taxa_names(data_fungi_mini)[1:20]
  sub <- prune_taxa(taxa_keep, data_fungi_mini)
  res <- compare_refseq_pq(data_fungi_mini, sub, verbose = FALSE)

  expect_equal(length(res$shared_names), ntaxa(sub))
  expect_gt(length(res$unique_names_1), 0)
  expect_equal(length(res$unique_names_2), 0)
  expect_false(is.na(res$mean_nn_kmer_dist_1))
  expect_true(is.na(res$mean_nn_kmer_dist_2))
  expect_gt(res$mean_nn_kmer_dist_1, 0)
})

test_that("compare_refseq_pq errors without refseq slot", {
  pq_no_refseq <- data_fungi_mini
  pq_no_refseq@refseq <- NULL

  expect_error(
    compare_refseq_pq(pq_no_refseq, data_fungi_mini, verbose = FALSE),
    "has no refseq slot"
  )
  expect_error(
    compare_refseq_pq(data_fungi_mini, pq_no_refseq, verbose = FALSE),
    "has no refseq slot"
  )
})

test_that("print method shows distance info when available", {
  taxa_keep <- taxa_names(data_fungi_mini)[1:20]
  sub <- prune_taxa(taxa_keep, data_fungi_mini)
  res <- compare_refseq_pq(data_fungi_mini, sub, verbose = FALSE)
  expect_output(print(res), "nearest-neighbor k-mer distance")
})

test_that("print method works for identical objects", {
  res <- compare_refseq_pq(data_fungi_mini, data_fungi_mini, verbose = FALSE)
  expect_output(print(res), "Reference Sequence Comparison")
})

test_that("compare_refseq_pq accepts list_phyloseq input", {
  sub <- prune_taxa(taxa_names(data_fungi_mini)[1:20], data_fungi_mini)
  lpq <- list_phyloseq(list(full = data_fungi_mini, subset = sub))
  res <- compare_refseq_pq(lpq, verbose = FALSE)

  expect_s3_class(res, "compare_refseq")
  expect_equal(res$name1, "full")
  expect_equal(res$name2, "subset")
  expect_equal(length(res$shared_names), ntaxa(sub))
})

test_that("compare_refseq_pq warns with >2 list_phyloseq objects", {
  lpq <- list_phyloseq(list(
    a = data_fungi_mini,
    b = data_fungi_mini,
    c = data_fungi_mini
  ))
  expect_message(
    compare_refseq_pq(lpq, verbose = FALSE),
    "Using the first two"
  )
})

test_that("compare_refseq_pq errors with <2 list_phyloseq objects", {
  lpq <- list_phyloseq(list(a = data_fungi_mini))
  expect_error(
    compare_refseq_pq(lpq, verbose = FALSE),
    "at least 2"
  )
})

test_that("compare_refseq_pq errors when physeq2 missing", {
  expect_error(
    compare_refseq_pq(data_fungi_mini, verbose = FALSE),
    "physeq2.*required"
  )
})

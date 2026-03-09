test_that("merge_lpq merges by names", {
  lpq <- list_phyloseq(list(
    fungi = data_fungi_mini,
    fungi2 = data_fungi_mini
  ))
  merged <- merge_lpq(lpq, match_by = "names", verbose = FALSE)
  expect_s4_class(merged, "phyloseq")
  expect_equal(phyloseq::nsamples(merged), 2L)
  expect_equal(
    as.character(phyloseq::sample_data(merged)$source_name),
    c("fungi", "fungi2")
  )
  expect_equal(phyloseq::ntaxa(merged), phyloseq::ntaxa(data_fungi_mini))
})

test_that("merge_lpq merges by refseq", {
  skip_if(
    is.null(phyloseq::refseq(data_fungi_mini, errorIfNULL = FALSE)),
    "data_fungi_mini has no refseq slot"
  )
  lpq <- list_phyloseq(list(
    fungi = data_fungi_mini,
    fungi2 = data_fungi_mini
  ))
  merged <- merge_lpq(lpq, match_by = "refseq", verbose = FALSE)
  expect_s4_class(merged, "phyloseq")
  expect_equal(phyloseq::nsamples(merged), 2L)
  expect_false(is.null(phyloseq::refseq(merged, errorIfNULL = FALSE)))
})

test_that("merge_lpq errors when refseq missing and match_by = 'refseq'", {
  pq_no_refseq <- phyloseq::phyloseq(
    phyloseq::otu_table(data_fungi_mini),
    phyloseq::tax_table(data_fungi_mini),
    phyloseq::sample_data(data_fungi_mini)
  )
  lpq <- list_phyloseq(list(a = pq_no_refseq, b = pq_no_refseq))
  expect_error(
    merge_lpq(lpq, match_by = "refseq", verbose = FALSE),
    "refseq slot"
  )
})

test_that("merge_lpq preserves total abundance per phyloseq", {
  lpq <- list_phyloseq(list(
    fungi = data_fungi_mini,
    fungi2 = data_fungi_mini
  ))
  merged <- merge_lpq(lpq, match_by = "names", verbose = FALSE)
  otu <- phyloseq::otu_table(merged)
  expected_sum <- sum(phyloseq::otu_table(data_fungi_mini))
  expect_equal(sum(otu[, "fungi"]), expected_sum)
  expect_equal(sum(otu[, "fungi2"]), expected_sum)
})

test_that("merge_lpq requires list_phyloseq input", {
  expect_error(
    merge_lpq(data_fungi_mini),
    "list_phyloseq"
  )
})

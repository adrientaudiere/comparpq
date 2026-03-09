# Tests for fake_creation.R
# Tests for functions that add fake taxa to phyloseq objects

# ==============================================================================
# Tests for add_shuffle_seq_pq
# ==============================================================================

test_that("add_shuffle_seq_pq adds correct number of fake taxa with n_fake", {
  set.seed(123)
  n_original <- ntaxa(data_fungi_mini)

  result <- add_shuffle_seq_pq(data_fungi_mini, n_fake = 10)

  expect_equal(ntaxa(result), n_original + 10)
})

test_that("add_shuffle_seq_pq adds correct proportion of fake taxa with prop_fake", {
  set.seed(123)
  n_original <- ntaxa(data_fungi_mini)
  prop <- 0.1

  result <- add_shuffle_seq_pq(data_fungi_mini, prop_fake = prop)

  expected_fake <- floor(prop * n_original)
  expect_equal(ntaxa(result), n_original + expected_fake)
})

test_that("add_shuffle_seq_pq fake taxa have NA in tax_table", {
  set.seed(123)
  result <- add_shuffle_seq_pq(data_fungi_mini, n_fake = 5)

  # Get fake taxa names
  fake_names <- grep("^fake_", taxa_names(result), value = TRUE)

  expect_equal(length(fake_names), 5)

  # Check that fake taxa have NA in tax_table
  fake_tax <- result@tax_table[fake_names, ]
  expect_true(all(is.na(fake_tax)))
})

test_that("add_shuffle_seq_pq fake taxa have 0 in otu_table", {
  set.seed(123)
  result <- add_shuffle_seq_pq(data_fungi_mini, n_fake = 5)

  # Get fake taxa names
  fake_names <- grep("^fake_", taxa_names(result), value = TRUE)

  # Check that fake taxa have 0 abundance
  if (taxa_are_rows(result)) {
    fake_otu <- result@otu_table[fake_names, ]
  } else {
    fake_otu <- result@otu_table[, fake_names]
  }

  expect_true(all(fake_otu == 0))
})

test_that("add_shuffle_seq_pq fails when neither n_fake nor prop_fake specified", {
  expect_error(
    add_shuffle_seq_pq(data_fungi_mini),
    "must specify either n_fake or prop_fake"
  )
})

test_that("add_shuffle_seq_pq fails when both n_fake and prop_fake specified", {
  expect_error(
    add_shuffle_seq_pq(data_fungi_mini, n_fake = 10, prop_fake = 0.1),
    "must specify either n_fake or prop_fake param, not both"
  )
})

test_that("add_shuffle_seq_pq respects custom prefix", {
  set.seed(123)
  result <- add_shuffle_seq_pq(
    data_fungi_mini,
    n_fake = 5,
    prefix = "shuffled_"
  )

  shuffled_names <- grep("^shuffled_", taxa_names(result), value = TRUE)

  expect_equal(length(shuffled_names), 5)
})

test_that("add_shuffle_seq_pq preserves original taxa", {
  set.seed(123)
  original_taxa <- taxa_names(data_fungi_mini)

  result <- add_shuffle_seq_pq(data_fungi_mini, n_fake = 5)

  expect_true(all(original_taxa %in% taxa_names(result)))
})

test_that("add_shuffle_seq_pq creates valid refseq for fake taxa", {
  set.seed(123)
  result <- add_shuffle_seq_pq(data_fungi_mini, n_fake = 5)

  # Get fake taxa names
  fake_names <- grep("^fake_", taxa_names(result), value = TRUE)

  # Check refseq exists for fake taxa
  fake_seqs <- result@refseq[fake_names]

  expect_equal(length(fake_seqs), 5)
  expect_s4_class(fake_seqs, "DNAStringSet")
})

# ==============================================================================
# Tests for add_external_seq_pq
# ==============================================================================

test_that("add_external_seq_pq adds external sequences", {
  ext_seqs <- Biostrings::readDNAStringSet(
    system.file("extdata/ex_little.fasta", package = "MiscMetabar")
  )

  n_original <- ntaxa(data_fungi_mini)

  result <- add_external_seq_pq(data_fungi_mini, ext_seqs)

  expect_equal(ntaxa(result), n_original + length(ext_seqs))
})

test_that("add_external_seq_pq external taxa have NA in tax_table", {
  ext_seqs <- Biostrings::readDNAStringSet(
    system.file("extdata/ex_little.fasta", package = "MiscMetabar")
  )

  result <- add_external_seq_pq(data_fungi_mini, ext_seqs)

  # Get external taxa names
  ext_names <- grep("^external_", taxa_names(result), value = TRUE)

  expect_equal(length(ext_names), length(ext_seqs))

  # Check that external taxa have NA in tax_table
  ext_tax <- result@tax_table[ext_names, ]
  expect_true(all(is.na(ext_tax)))
})

test_that("add_external_seq_pq external taxa have 0 in otu_table", {
  ext_seqs <- Biostrings::readDNAStringSet(
    system.file("extdata/ex_little.fasta", package = "MiscMetabar")
  )

  result <- add_external_seq_pq(data_fungi_mini, ext_seqs)

  # Get external taxa names
  ext_names <- grep("^external_", taxa_names(result), value = TRUE)

  # Check that external taxa have 0 abundance
  if (taxa_are_rows(result)) {
    ext_otu <- result@otu_table[ext_names, ]
  } else {
    ext_otu <- result@otu_table[, ext_names]
  }

  expect_true(all(ext_otu == 0))
})

test_that("add_external_seq_pq respects custom prefix", {
  ext_seqs <- Biostrings::readDNAStringSet(
    system.file("extdata/ex_little.fasta", package = "MiscMetabar")
  )

  result <- add_external_seq_pq(data_fungi_mini, ext_seqs, prefix = "added_")

  added_names <- grep("^added_", taxa_names(result), value = TRUE)

  expect_equal(length(added_names), length(ext_seqs))
})

test_that("add_external_seq_pq preserves original taxa", {
  ext_seqs <- Biostrings::readDNAStringSet(
    system.file("extdata/ex_little.fasta", package = "MiscMetabar")
  )

  original_taxa <- taxa_names(data_fungi_mini)

  result <- add_external_seq_pq(data_fungi_mini, ext_seqs)

  expect_true(all(original_taxa %in% taxa_names(result)))
})

test_that("add_external_seq_pq preserves sample_data structure", {
  ext_seqs <- Biostrings::readDNAStringSet(
    system.file("extdata/ex_little.fasta", package = "MiscMetabar")
  )

  result <- add_external_seq_pq(data_fungi_mini, ext_seqs)

  expect_equal(nsamples(result), nsamples(data_fungi_mini))
  expect_equal(sample_names(result), sample_names(data_fungi_mini))
})

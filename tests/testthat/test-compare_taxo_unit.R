# Tests for compare_taxo_unit.R: scoring against a per-unit truth, where every
# scored unit falls in exactly one cell of the confusion matrix.

# A hand-made object: 4 real units, 2 shuffled controls, 2 external controls.
# The tax_table holds three assignations of the Genus rank and one of Kingdom.
unit_mock <- function() {
  taxa <- c(
    "u_correct",
    "u_wrong",
    "u_na",
    "u_no_truth",
    "u_tie",
    "fake_1",
    "fake_2",
    "external_1",
    "external_2"
  )
  otu <- phyloseq::otu_table(
    matrix(1, nrow = 1, ncol = length(taxa), dimnames = list("s1", taxa)),
    taxa_are_rows = FALSE
  )
  tax <- phyloseq::tax_table(as.matrix(data.frame(
    row.names = taxa,
    Kingdom_good = c(
      "Fungi",
      "Fungi",
      NA,
      "Fungi",
      "Fungi",
      NA,
      "Fungi",
      "Fungi",
      "Viridiplantae"
    ),
    Genus_good = c(
      "Alpha",
      "Beta",
      NA,
      "Delta",
      "Gamma",
      NA,
      "Zeta",
      NA,
      "Eta"
    ),
    Genus_accepted_name = c(
      "Alpha_new",
      "Beta",
      NA,
      NA,
      "Gamma",
      NA,
      NA,
      NA,
      NA
    ),
    Genus_silent = rep(NA_character_, length(taxa))
  )))
  phyloseq::phyloseq(otu, tax)
}

unit_truth <- function() {
  data.frame(
    taxon = c("u_correct", "u_wrong", "u_na", "u_tie"),
    Kingdom = "Fungi",
    Genus = c("Alpha", "Alpha", "Alpha", "Gamma"),
    Genus_accepted = c("Alpha_new", "Alpha_new", "Alpha_new", NA),
    Species = c("Alpha_one", "Alpha_one", "Alpha_one", NA),
    truth_depth = c("Species", "Species", "Species", "Genus")
  )
}

ranks <- c("Kingdom", "Genus", "Species")

test_that("every scored unit falls in exactly one cell", {
  res <- tc_metrics_unit_vec(
    unit_mock(),
    taxonomic_rank = "Genus_good",
    truth = unit_truth(),
    rank = "Genus",
    truth_ranks = ranks,
    external_scoring = "matrix",
    verbose = FALSE
  )
  # 5 real + 2 shuffled + 2 external = 9 units, none left out at Genus.
  expect_equal(res$TP + res$FP + res$FN + res$TN, 9)
  expect_equal(res$TP, 2) # u_correct (Alpha) and u_tie (Gamma)
  expect_equal(res$FN, 1) # u_na
  # u_wrong (Beta), u_no_truth (Delta), fake_2 (Zeta), external_2 (Eta)
  expect_equal(res$FP, 4)
  expect_equal(res$TN, 2) # fake_1 and external_1, left NA
  expect_equal(res$TNR, 2 / 6)
  expect_equal(res$NA_fake, 0.5)
  expect_equal(res$ctrl_assigned_fake, 1)
  expect_equal(res$ctrl_assigned_ext, 1)
  expect_equal(res$n_real, 5)
  expect_equal(res$n_left_out, 0)
  expect_equal(res$misassign_seq, 2 / 5) # u_wrong and u_no_truth
  expect_equal(res$NA_real, 1 / 5)
})

test_that("the accepted name counts as correct too", {
  res <- tc_metrics_unit_vec(
    unit_mock(),
    taxonomic_rank = "Genus_accepted_name",
    truth = unit_truth(),
    rank = "Genus",
    truth_ranks = ranks,
    external_scoring = "matrix",
    verbose = FALSE
  )
  # u_correct is named with the accepted name, u_wrong keeps Beta.
  expect_equal(res$TP, 2)
  expect_equal(res$FP, 1)
  expect_equal(res$TN, 4)
})

test_that("a unit is left out of the matrix below its truth", {
  res <- tc_metrics_unit_vec(
    unit_mock(),
    taxonomic_rank = "Genus_good",
    truth = unit_truth(),
    rank = "Species",
    truth_ranks = ranks,
    external_scoring = "matrix",
    verbose = FALSE
  )
  # u_tie stops at Genus: it is not scored at Species.
  expect_equal(res$n_left_out, 1)
  expect_equal(res$n_real, 4)
  expect_equal(res$TP + res$FP + res$FN + res$TN, 8)
})

test_that("external controls stay aside when the database can hold them", {
  aside <- tc_metrics_unit_vec(
    unit_mock(),
    taxonomic_rank = "Genus_good",
    truth = unit_truth(),
    rank = "Genus",
    truth_ranks = ranks,
    external_scoring = "aside",
    verbose = FALSE
  )
  expect_equal(aside$TP + aside$FP + aside$FN + aside$TN, 7)
  expect_equal(aside$TN, 1) # fake_1 only
  expect_equal(aside$FP, 3) # external_2 no longer counted
  expect_equal(aside$ctrl_assigned_ext, 1) # still reported
})

test_that("a column assigning nothing gives TN on the controls and FN on the real units", {
  res <- tc_metrics_unit_vec(
    unit_mock(),
    taxonomic_rank = "Genus_silent",
    truth = unit_truth(),
    rank = "Genus",
    truth_ranks = ranks,
    external_scoring = "matrix",
    verbose = FALSE
  )
  expect_equal(res$TP, 0)
  expect_equal(res$FP, 0)
  expect_equal(res$FN, 5)
  expect_equal(res$TN, 4)
  expect_equal(res$NA_real, 1)
  expect_equal(res$NA_fake, 1)
  expect_equal(res$MCC, 0) # zero margin, Chicco & Jurman 2020
  expect_equal(res$TNR, 1)
})

test_that("tc_metrics_unit adds ext_fungi and ext_correct", {
  ranks_df <- data.frame(good = c("Kingdom_good", "Genus_good", "Genus_good"))
  external_truth <- data.frame(
    taxon = c("external_1", "external_2"),
    Kingdom = c("Viridiplantae", "Viridiplantae"),
    Genus = c("Pinus", "Eta"),
    Species = c("Pinus_pinea", "Eta_one")
  )
  res <- tc_metrics_unit(
    unit_mock(),
    ranks_df = ranks_df,
    truth = unit_truth(),
    truth_ranks = ranks,
    external_truth = external_truth,
    external_scoring = "aside"
  )
  expect_setequal(
    colnames(res),
    c("method_db", "tax_level", "metrics", "values")
  )
  expect_setequal(unique(res$tax_level), ranks)
  ext_fungi <- res$values[res$metrics == "ext_fungi"]
  expect_length(ext_fungi, 1) # kingdom rank only
  expect_equal(ext_fungi, 0.5) # external_1 called Fungi
  ext_correct <- res$values[
    res$metrics == "ext_correct" & res$tax_level == "Genus"
  ]
  expect_equal(ext_correct, 0.5) # external_2 named Eta
  expect_equal(
    res$values[res$metrics == "ext_correct" & res$tax_level == "Kingdom"],
    0.5
  ) # external_2 named Viridiplantae
})

test_that("tc_metrics_unit counts ext_correct over the controls with a truth", {
  ranks_df <- data.frame(good = c("Kingdom_good", "Genus_good", "Genus_good"))
  external_truth <- data.frame(
    taxon = c("external_1", "external_2"),
    Kingdom = "Viridiplantae",
    Genus = c(NA, "Eta"),
    Species = c(NA, "Eta_one")
  )
  res <- tc_metrics_unit(
    unit_mock(),
    ranks_df = ranks_df,
    truth = unit_truth(),
    truth_ranks = ranks,
    external_truth = external_truth,
    external_scoring = "aside"
  )
  # external_1 has no genus truth: left out, external_2 named Eta.
  expect_equal(
    res$values[res$metrics == "ext_correct" & res$tax_level == "Genus"],
    1
  )
})

test_that("tc_metrics_unit can restrict ext_correct to the kingdom", {
  ranks_df <- data.frame(good = c("Kingdom_good", "Genus_good", "Genus_good"))
  external_truth <- data.frame(
    taxon = c("external_1", "external_2"),
    Kingdom = "Viridiplantae",
    Genus = c("Pinus", "Eta"),
    Species = c("Pinus_pinea", "Eta_one")
  )
  res <- tc_metrics_unit(
    unit_mock(),
    ranks_df = ranks_df,
    truth = unit_truth(),
    truth_ranks = ranks,
    external_truth = external_truth,
    external_correct_ranks = "Kingdom",
    external_scoring = "aside"
  )
  expect_equal(sum(res$metrics == "ext_correct"), 1)
  expect_equal(unique(res$tax_level[res$metrics == "ext_correct"]), "Kingdom")
})

test_that("tc_metrics_unit refuses a ranks_df that does not match the ranks", {
  expect_error(
    tc_metrics_unit(
      unit_mock(),
      ranks_df = data.frame(good = c("Kingdom_good", "Genus_good")),
      truth = unit_truth(),
      truth_ranks = ranks
    ),
    "one row per rank"
  )
  expect_error(
    tc_metrics_unit_vec(
      unit_mock(),
      taxonomic_rank = "Genus_good",
      truth = unit_truth(),
      rank = "Phylum",
      truth_ranks = ranks
    ),
    "not one of"
  )
})

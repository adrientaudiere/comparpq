test_that("extract_species_epithet handles the documented formats", {
  expect_identical(extract_species_epithet("ilex", "Quercus"), "ilex")
  expect_identical(extract_species_epithet("Quercus ilex", "Quercus"), "ilex")
  expect_identical(extract_species_epithet("Quercus_ilex", "Quercus"), "ilex")
  expect_identical(
    extract_species_epithet("Quercus_ilex_var._rotundifolia", "Quercus"),
    "ilex"
  )
  expect_identical(extract_species_epithet("Quercus", "Quercus"), NA_character_)
  expect_identical(
    extract_species_epithet(NA_character_, "Quercus"),
    NA_character_
  )
  expect_identical(extract_species_epithet("", "Quercus"), NA_character_)
})

test_that("harmonize_sp_names_pq (verify = FALSE) keeps only the epithet", {
  pq <- suppressMessages(
    MiscMetabar::subset_taxa_pq(Glom_otu, phyloseq::taxa_sums(Glom_otu) > 8000)
  )
  mat <- as(pq@tax_table, "matrix")
  n <- nrow(mat)
  mat <- cbind(
    mat,
    Genus_EUK = rep("Quercus", n),
    Species_EUK = rep(
      c("Quercus_ilex", "Quercus robur", "suber"),
      length.out = n
    )
  )
  pq@tax_table <- phyloseq::tax_table(mat)

  out <- harmonize_sp_names_pq(pq, suffixes = "_EUK", verify = FALSE)
  sp <- as.character(out@tax_table[, "Species_EUK"])
  expect_true(all(sp %in% c("ilex", "robur", "suber")))
})

test_that("harmonize_sp_names_pq normalises the suffix (no leading underscore)", {
  pq <- suppressMessages(
    MiscMetabar::subset_taxa_pq(Glom_otu, phyloseq::taxa_sums(Glom_otu) > 8000)
  )
  mat <- as(pq@tax_table, "matrix")
  n <- nrow(mat)
  mat <- cbind(
    mat,
    Genus_EUK = rep("Quercus", n),
    Species_EUK = rep("Quercus_ilex", n)
  )
  pq@tax_table <- phyloseq::tax_table(mat)

  out <- harmonize_sp_names_pq(pq, suffixes = "EUK", verify = FALSE)
  expect_true(all(as.character(out@tax_table[, "Species_EUK"]) == "ilex"))
})

test_that("harmonize_sp_names_pq aborts on verify = TRUE without taxinfo", {
  skip_if(requireNamespace("taxinfo", quietly = TRUE))
  pq <- suppressMessages(
    MiscMetabar::subset_taxa_pq(Glom_otu, phyloseq::taxa_sums(Glom_otu) > 8000)
  )
  expect_error(
    harmonize_sp_names_pq(pq, suffixes = "_EUK", verify = TRUE),
    "taxinfo"
  )
})

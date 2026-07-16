test_that("count_taxo_congruence returns a 5-category tibble", {
  res <- count_taxo_congruence(Glom_otu, "Genus", "Genus__eukaryome_Glomero")

  expect_s3_class(res, "tbl_df")
  expect_identical(nrow(res), 5L)
  expect_identical(
    colnames(res),
    c("category", "n_asv", "n_seq", "pct_asv", "pct_seq")
  )
  expect_true("both_equal" %in% res$category)
  expect_true("only_Genus" %in% res$category)
  expect_true("only_Genus__eukaryome_Glomero" %in% res$category)
})

test_that("count_taxo_congruence counts and percentages are coherent", {
  res <- count_taxo_congruence(Glom_otu, "Genus", "Genus__eukaryome_Glomero")

  expect_equal(sum(res$n_asv), phyloseq::ntaxa(Glom_otu))
  expect_equal(sum(res$pct_asv), 100, tolerance = 0.05)
  expect_equal(sum(res$pct_seq), 100, tolerance = 0.05)
  expect_true(all(res$n_asv >= 0))
})

test_that("count_taxo_congruence treats empty and 'NA' strings as missing", {
  pq <- Glom_otu
  mat <- as(pq@tax_table, "matrix")
  mat[1, "Genus"] <- ""
  mat[2, "Genus"] <- "NA"
  pq@tax_table <- phyloseq::tax_table(mat)

  res <- count_taxo_congruence(pq, "Genus", "Genus__eukaryome_Glomero")
  only_g2 <- res$n_asv[res$category == "only_Genus__eukaryome_Glomero"]
  expect_true(only_g2 >= 2)
})

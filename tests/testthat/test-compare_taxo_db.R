make_small_pq <- function() {
  suppressMessages(
    MiscMetabar::subset_taxa_pq(Glom_otu, phyloseq::taxa_sums(Glom_otu) > 8000)
  )
}

test_that("compare_taxo_db returns the expected bundle", {
  pq <- make_small_pq()
  res <- suppressWarnings(
    compare_taxo_db(pq, "Genus", "Genus__eukaryome_Glomero")
  )

  expect_type(res, "list")
  expect_identical(
    names(res),
    c("congruence", "contingency", "bar", "sankey", "circle", "heatmap")
  )
  expect_s3_class(res$bar, "ggplot")
  expect_s3_class(res$sankey, "ggplot")
  expect_s3_class(res$heatmap, "ggplot")
  expect_type(res$circle, "closure")
  expect_true(is.data.frame(res$contingency))
})

test_that("compare_taxo_db contingency respects min_n and ordering", {
  pq <- make_small_pq()
  res <- suppressWarnings(
    compare_taxo_db(pq, "Genus", "Genus__eukaryome_Glomero", min_n = 1)
  )
  expect_true(all(res$contingency$n > 1))
  expect_false(is.unsorted(rev(res$contingency$n)))
})

test_that("compare_taxo_db resolves default rank_3 from column naming", {
  pq <- make_small_pq()
  res <- suppressWarnings(
    compare_taxo_db(
      pq,
      "Genus__eukaryome_Glomero",
      "Genus__eukaryome_Glomero"
    )
  )
  expect_s3_class(res$bar, "ggplot")
})

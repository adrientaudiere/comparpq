test_that("plot_congruence_counts works with explicit ranks_1/ranks_2", {
  p <- plot_congruence_counts(
    Glom_otu,
    ranks_1 = c("Order", "Family", "Genus"),
    ranks_2 = c(
      "Order__eukaryome_Glomero",
      "Family__eukaryome_Glomero",
      "Genus__eukaryome_Glomero"
    )
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_congruence_counts works in suffix mode", {
  pq <- mock_two_db_pq()
  p <- plot_congruence_counts(pq, suffix_1 = "_A", suffix_2 = "_B")
  expect_s3_class(p, "ggplot")
})

test_that("plot_congruence_counts honours n_seq", {
  p <- plot_congruence_counts(
    Glom_otu,
    ranks_1 = "Genus",
    ranks_2 = "Genus__eukaryome_Glomero",
    n_seq = TRUE
  )
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$y, "% of sequences")
})

test_that("plot_congruence_counts errors on mismatched lengths", {
  expect_error(
    plot_congruence_counts(
      Glom_otu,
      ranks_1 = c("Genus", "Family"),
      ranks_2 = "Genus__eukaryome_Glomero"
    ),
    "same length"
  )
})

test_that("plot_congruence_counts errors when only one suffix given", {
  expect_error(
    plot_congruence_counts(Glom_otu, suffix_1 = "_EUK"),
    "provided together"
  )
})

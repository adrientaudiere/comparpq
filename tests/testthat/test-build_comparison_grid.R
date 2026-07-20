test_that("build_comparison_grid enumerates rank x database pairs", {
  pq <- mock_two_db_pq()
  grid <- build_comparison_grid(list(m1 = pq))

  expect_s3_class(grid, "tbl_df")
  expect_identical(
    colnames(grid),
    c("physeq", "rank", "db_1", "col_1", "db_2", "col_2")
  )
  # 3 default ranks x 3 choose 2 database pairs (A, B, _eukaryome_Glomero)
  expect_equal(nrow(grid), 9)
  expect_true(all(grid$physeq == "m1"))
  expect_setequal(unique(grid$rank), c("Order", "Family", "Genus"))
  expect_true(all(grid$db_1 != grid$db_2))
  expect_true(all(grid$col_1 %in% colnames(pq@tax_table)))
  expect_true(all(grid$col_2 %in% colnames(pq@tax_table)))
})

test_that("build_comparison_grid respects the ranks argument", {
  pq <- mock_two_db_pq()
  grid <- build_comparison_grid(list(m1 = pq), ranks = "Species")

  expect_equal(nrow(grid), 3)
  expect_true(all(grid$rank == "Species"))
})

test_that("build_comparison_grid handles several phyloseq objects", {
  pq <- mock_two_db_pq()
  grid <- build_comparison_grid(list(first = pq, second = pq))

  expect_equal(nrow(grid), 18)
  expect_setequal(unique(grid$physeq), c("first", "second"))
})

test_that("build_comparison_grid auto-names unnamed lists", {
  pq <- mock_two_db_pq()
  expect_message(
    grid <- build_comparison_grid(list(pq)),
    "auto-named"
  )
  expect_true(all(grid$physeq == "pq1"))
})

test_that("build_comparison_grid accepts a list_phyloseq object", {
  pq <- mock_two_db_pq()
  lpq <- suppressMessages(
    list_phyloseq(list(m1 = pq), compute_dist = FALSE, verbose = FALSE)
  )
  grid <- build_comparison_grid(lpq)

  expect_equal(nrow(grid), 9)
  expect_true(all(grid$physeq == "m1"))
})

test_that("build_comparison_grid returns an empty grid without pairs", {
  pq <- mock_two_db_pq()
  mat <- as(pq@tax_table, "matrix")
  pq@tax_table <- phyloseq::tax_table(mat[, c("Order", "Family", "Genus")])

  grid <- build_comparison_grid(list(m1 = pq))
  expect_s3_class(grid, "tbl_df")
  expect_equal(nrow(grid), 0)
  expect_identical(
    colnames(grid),
    c("physeq", "rank", "db_1", "col_1", "db_2", "col_2")
  )
})

test_that("build_comparison_grid validates its input", {
  pq <- mock_two_db_pq()
  expect_error(build_comparison_grid(list()), "non-empty")
  expect_error(
    build_comparison_grid(list(m1 = pq, bad = data.frame())),
    "must be phyloseq objects"
  )
})

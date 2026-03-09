# Tests for gg_hill_lpq.R

create_test_lpq <- function() {
  list_phyloseq(
    list(run1 = data_fungi, run2 = data_fungi),
    same_bioinfo_pipeline = FALSE
  )
}

test_that("gg_hill_lpq returns a ggplot", {
  lpq <- create_test_lpq()
  p <- gg_hill_lpq(lpq, verbose = FALSE)
  expect_s3_class(p, "gg")
})

test_that("gg_hill_lpq works with a single hill scale", {
  lpq <- create_test_lpq()
  p <- gg_hill_lpq(lpq, hill_scales = 0, verbose = FALSE)
  expect_s3_class(p, "gg")
})

test_that("gg_hill_lpq works with three phyloseq objects", {
  lpq <- list_phyloseq(
    list(run1 = data_fungi, run2 = data_fungi, run3 = data_fungi),
    same_bioinfo_pipeline = FALSE
  )
  p <- gg_hill_lpq(lpq, verbose = FALSE)
  expect_s3_class(p, "gg")
})

test_that("gg_hill_lpq respects the pairs argument", {
  lpq <- list_phyloseq(
    list(a = data_fungi, b = data_fungi, c = data_fungi),
    same_bioinfo_pipeline = FALSE
  )
  p <- gg_hill_lpq(lpq, pairs = list(c(1L, 2L)), verbose = FALSE)
  expect_s3_class(p, "gg")
})

test_that("gg_hill_lpq errors with fewer than 2 phyloseq objects", {
  lpq <- list_phyloseq(list(only = data_fungi))
  expect_error(gg_hill_lpq(lpq), "at least 2")
})

test_that("gg_hill_lpq errors with invalid cor_method", {
  lpq <- create_test_lpq()
  expect_error(gg_hill_lpq(lpq, cor_method = "bad"), "cor_method must be")
})

test_that("gg_hill_lpq works with spearman correlation", {
  lpq <- create_test_lpq()
  p <- gg_hill_lpq(lpq, cor_method = "spearman", verbose = FALSE)
  expect_s3_class(p, "gg")
})

test_that("gg_hill_lpq works without smooth and without 1to1 line", {
  lpq <- create_test_lpq()
  p <- gg_hill_lpq(lpq, add_smooth = FALSE, add_1to1 = FALSE, verbose = FALSE)
  expect_s3_class(p, "gg")
})

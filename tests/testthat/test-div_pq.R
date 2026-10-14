skip_on_cran()

test_that("div_pq() accepts fact and deprecates modality", {
  pq_sub <- prune_samples(
    sample_names(data_fungi_mini)[1:20],
    data_fungi_mini
  )
  res_fact <- div_pq(
    pq_sub,
    fact = "Height",
    indices = c("shannon", "simpson"),
    q = c(0, 1)
  )
  lifecycle::expect_deprecated(
    res_mod <- div_pq(
      pq_sub,
      modality = "Height",
      indices = c("shannon", "simpson"),
      q = c(0, 1)
    )
  )
  expect_identical(res_fact, res_mod)
  expect_true("Height" %in% colnames(res_fact))
})

test_that("div_pq() aggregation and significance work with fact", {
  skip_if_not_installed("multcompView")
  pq_sub <- prune_samples(
    sample_names(data_fungi_mini)[1:20],
    data_fungi_mini
  )
  res_agg <- div_pq(
    pq_sub,
    fact = "Height",
    indices = "shannon",
    aggregate = TRUE,
    significance = TRUE
  )
  expect_true("Height" %in% colnames(res_agg))
  expect_true("shannon_(kruskal.p)" %in% colnames(res_agg))
  expect_true("shannon_(tukey.letters)" %in% colnames(res_agg))
})

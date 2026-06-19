# Extracted from test-estim_pq.R:434

# prequel ----------------------------------------------------------------------
library("divent")
create_test_pq <- function() {
  pq <- subset_samples(data_fungi, !is.na(Height))
  clean_pq(pq)
}
create_test_lpq <- function() {
  pq <- create_test_pq()
  list_phyloseq(
    list(original = pq, copy = pq),
    same_bioinfo_pipeline = FALSE
  )
}

# test -------------------------------------------------------------------------
skip_if_not_installed("dabestr")
lpq <- create_test_lpq()
res <- estim_diff_lpq(
  lpq,
  fact = "Height",
  q = c(0),
  resamples = 100,
  verbose = FALSE
)
expect_output(print(res), "list_phyloseq")

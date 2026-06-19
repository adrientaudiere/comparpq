# Extracted from test-estim_pq.R:501

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
pq <- create_test_pq()
sam <- sample_data(pq)
sam$lib_size <- sample_sums(pq)
sample_data(pq) <- sam
lpq <- list_phyloseq(
  list(original = pq, copy = pq),
  same_bioinfo_pipeline = FALSE
)
set.seed(42)
res <- estim_cor_lpq(
  lpq,
  variable = "lib_size",
  q = c(0),
  resamples = 100,
  verbose = FALSE
)
expect_output(print(res), "list_phyloseq")

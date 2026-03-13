# Estimation statistics for numeric correlation on a list_phyloseq

Applies
[`estim_cor_pq()`](https://adrientaudiere.github.io/comparpq/reference/estim_cor_pq.md)
to each phyloseq object in a list_phyloseq and combines the results.

## Usage

``` r
estim_cor_lpq(x, variable, ..., verbose = TRUE)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- variable:

  (character, required) The name of a numeric column in `sample_data`.
  Must be present in **all** phyloseq objects.

- ...:

  Additional arguments passed to
  [`estim_cor_pq()`](https://adrientaudiere.github.io/comparpq/reference/estim_cor_pq.md).

- verbose:

  (logical, default TRUE) If TRUE, print progress messages.

## Value

A list of class `"estim_cor_lpq_result"` with components:

- results:

  A named list of `estim_cor_pq_result` objects (one per phyloseq)

- correlations:

  A tibble combining all correlations with an additional `name` column

- regressions:

  A tibble combining all regressions with an additional `name` column

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## See also

[`estim_cor_pq()`](https://adrientaudiere.github.io/comparpq/reference/estim_cor_pq.md),
[`estim_diff_lpq()`](https://adrientaudiere.github.io/comparpq/reference/estim_diff_lpq.md)

## Examples

``` r
if (FALSE) { # \dontrun{
lpq <- list_phyloseq(
  list(
    fungi = data_fungi,
    fungi_clust = postcluster_pq(data_fungi)
  ),
  same_bioinfo_pipeline = FALSE
)

# Assuming a numeric variable exists in sample_data
results <- estim_cor_lpq(lpq, variable = "lib_size")
results$correlations
} # }
```

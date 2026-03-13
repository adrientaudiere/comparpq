# Estimation statistics for numeric variable correlation on a phyloseq object

Computes diversity metrics (Hill numbers by default) per sample and
assesses their relationship with a numeric variable using bootstrap
confidence intervals for correlation coefficients and regression slopes.

## Usage

``` r
estim_cor_pq(
  physeq,
  variable,
  hill_scales = c(0, 1, 2),
  custom_fn = NULL,
  method = "pearson",
  resamples = 5000,
  ci = 95,
  na_remove = TRUE
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- variable:

  (character, required) The name of a numeric column in `sample_data`.

- hill_scales:

  (numeric vector, default `c(0, 1, 2)`) The q values for Hill number
  computation.

- custom_fn:

  (function, default NULL) An optional custom diversity function (see
  [`estim_diff_pq()`](https://adrientaudiere.github.io/comparpq/reference/estim_diff_pq.md)
  for details).

- method:

  (character, default `"pearson"`) Correlation method. One of
  `"pearson"`, `"spearman"`, `"kendall"`.

- resamples:

  (integer, default 5000) Number of bootstrap resamples.

- ci:

  (numeric, default 95) Confidence interval level (0-100).

- na_remove:

  (logical, default TRUE) If TRUE, samples with NA in `variable` are
  removed.

## Value

A list of class `"estim_cor_pq_result"` with components:

- data:

  The diversity data.frame used for analysis

- correlations:

  A tibble with columns: `metric`, `estimate`, `ci_lower`, `ci_upper`,
  `method`, `pvalue`

- regressions:

  A tibble with columns: `metric`, `intercept`, `slope`,
  `slope_ci_lower`, `slope_ci_upper`

- plots:

  A named list of ggplot2 scatter plots with regression line and
  bootstrap CI ribbon (one per metric)

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## See also

[`estim_diff_pq()`](https://adrientaudiere.github.io/comparpq/reference/estim_diff_pq.md),
[`estim_cor_lpq()`](https://adrientaudiere.github.io/comparpq/reference/estim_cor_lpq.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(phyloseq)
data("data_fungi", package = "MiscMetabar")

# Add a numeric variable for demonstration
sam <- sample_data(data_fungi)
sam$lib_size <- sample_sums(data_fungi)
sample_data(data_fungi) <- sam

res <- estim_cor_pq(data_fungi, variable = "lib_size")
res
res$plots$Hill_0
res$correlations
} # }
```

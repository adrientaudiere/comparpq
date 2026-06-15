# Estimation statistics for categorical comparisons on a phyloseq object

Computes diversity metrics (Hill numbers by default) per sample and
compares them across groups defined by a categorical variable using
estimation statistics (effect sizes + bootstrap confidence intervals)
via the dabestr package.

This approach replaces traditional p-value-based hypothesis testing with
Gardner-Altman or Cumming estimation plots, following the estimation
statistics framework (Ho et al. 2019).

## Usage

``` r
estim_diff_pq(
  physeq,
  fact,
  q = c(0, 1, 2),
  custom_fn = NULL,
  effect_type = "cohens_d",
  idx = NULL,
  ci = 95,
  resamples = 5000,
  na_remove = TRUE,
  ...
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- fact:

  (character, required) The name of a categorical column in
  `sample_data` to use as the grouping factor.

- q:

  (numeric vector, default `c(0, 1, 2)`) The q values for Hill number
  computation: 0 = richness, 1 = Shannon exponential, 2 = inverse
  Simpson.

- custom_fn:

  (function, default NULL) An optional custom diversity function. Must
  take a phyloseq object and return a named numeric vector (names =
  sample names) or a data.frame with one row per sample. If provided,
  `q` is ignored.

- effect_type:

  (character, default `"mean_diff"`) The type of effect size to compute.
  One of: `"mean_diff"`, `"median_diff"`, `"cohens_d"`, `"hedges_g"`,
  `"cliffs_delta"`.

- idx:

  (list or character vector, default NULL) The group ordering for
  comparisons. If NULL, uses factor levels with first level as control.
  For 2 groups: `c("Control", "Treatment")`. For 3+ groups:
  `list(c("Ctrl", "T1", "T2"))`.

- ci:

  (numeric, default 95) Confidence interval level (0-100).

- resamples:

  (integer, default 5000) Number of bootstrap resamples.

- na_remove:

  (logical, default TRUE) If TRUE, samples with NA in `fact` are removed
  before analysis.

- ...:

  Additional arguments passed to dabestr plotting functions.

## Value

A list of class `"estim_diff_pq_result"` with components:

- data:

  The diversity data.frame used for analysis

- dabest_objects:

  A named list of dabestr objects (one per metric)

- plots:

  A named list of dabestr plots (one per metric)

- summary:

  A tibble summarizing all effect sizes and CIs with columns: `metric`,
  `comparison`, `effect_size`, `ci_lower`, `ci_upper`,
  `pvalue_permtest`, `pvalue_welch`, `pvalue_mann_whitney`

- effect_type:

  The effect size type used

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

The function uses the dabestr package to produce estimation plots. For
two groups, Gardner-Altman plots are produced (`float_contrast = TRUE`).
For three or more groups, Cumming plots are used
(`float_contrast = FALSE`).

## References

Ho, J., Tumkaya, T., Aryal, S., Choi, H., & Claridge-Chang, A. (2019).
Moving beyond P values: data analysis with estimation graphics. *Nature
Methods*, 16(7), 565-566.

## See also

[`estim_cor_pq()`](https://adrientaudiere.github.io/comparpq/reference/estim_cor_pq.md),
[`estim_diff_lpq()`](https://adrientaudiere.github.io/comparpq/reference/estim_diff_lpq.md),
[`adonis_lpq()`](https://adrientaudiere.github.io/comparpq/reference/adonis_lpq.md)

## Examples

``` r
library(phyloseq)
data("data_fungi", package = "MiscMetabar")

pq <- subset_samples(data_fungi, !is.na(Height))
pq <- clean_pq(pq)
#> Cleaning suppress 144 taxa and 0 samples.

res <- estim_diff_pq(pq, fact = "Height")
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
res
#> $data
#>                                   Hill_0    Hill_1    Hill_2
#> A10-005-B_S188_MERGED.fastq.gz        76  8.001711  3.178277
#> A10-005-H_S189_MERGED.fastq.gz        97  2.708544  1.908804
#> A10-005-M_S190_MERGED.fastq.gz        63 18.977884 13.213389
#> A12-007-B_S2_MERGED.fastq.gz          56  9.043012  5.231148
#> AB29-ABMX-H_S6_MERGED.fastq.gz        80 10.257713  5.505835
#> AD26-005-B_S9_MERGED.fastq.gz        236  9.037148  3.358779
#> AD26-005-H_S10_MERGED.fastq.gz       125  4.949292  2.385584
#> AD26-005-M_S11_MERGED.fastq.gz       141 27.836759 16.469898
#> AD30-ABMX-M_S12_MERGED.fastq.gz      100  8.470705  4.483098
#> AD32-007-M_S13_MERGED.fastq.gz       118  8.435683  3.324192
#> ADABM30X-B_S14_MERGED.fastq.gz       123 15.365729 10.768405
#> ADABM30X-H_S15_MERGED.fastq.gz       159 18.060943 10.736997
#> ADABM30X-M_S16_MERGED.fastq.gz        81 11.528659  6.839672
#> B18-006-B_S19_MERGED.fastq.gz         35  4.664783  2.465271
#> BA17-050-B_S21_MERGED.fastq.gz       189 15.069287  7.464093
#> BB19-006-H_S22_MERGED.fastq.gz       149  8.263810  4.534396
#> BB6-019-B_S23_MERGED.fastq.gz        110  6.039948  3.163498
#> BB6-019-H_S24_MERGED.fastq.gz        141  7.174113  2.623084
#> BB6-019-M_S25_MERGED.fastq.gz        114 26.050805 15.082754
#> BE9-006-B_S27_MERGED.fastq.gz         20  5.983268  3.372952
#> BE9-006-H_S28_MERGED.fastq.gz         14  5.869285  3.773643
#> BE9-006-M_S29_MERGED.fastq.gz         27 12.262107  8.902889
#> BG7-010-B_S30_MERGED.fastq.gz         47 16.372452  6.686758
#> BG7-010-H_S31_MERGED.fastq.gz         96 11.189039  4.419449
#> BG7-010-M_S32_MERGED.fastq.gz         59 15.585219  6.620295
#> BJ17-007-M_S34_MERGED.fastq.gz       192 10.933159  5.789013
#> BL7-006-B_S36_MERGED.fastq.gz         91 20.436676  7.037519
#> BL7-006-H_S37_MERGED.fastq.gz         60 23.752281 11.122957
#> BL7-006-M_S38_MERGED.fastq.gz         42 17.788122  9.843910
#> BP11-001-B_S43_MERGED.fastq.gz       154  6.672521  3.173428
#> BP11-001-H_S44_MERGED.fastq.gz       126  8.080340  2.946738
#> BP11-001-M_S45_MERGED.fastq.gz        80  3.923054  1.957611
#> BP12-025-B_S46_MERGED.fastq.gz        62  7.801954  4.053323
#> BQ4-018-B_S49_MERGED.fastq.gz        239  8.071347  4.608665
#> BQ4-018-H_S50_MERGED.fastq.gz         84  9.209574  4.825529
#> BQ4-018-M_S51_MERGED.fastq.gz        206 16.210572  8.277339
#> BT-006-M_S55_MERGED.fastq.gz          68  9.743600  6.005182
#> BV11-002-B_S57_MERGED.fastq.gz       120  6.975027  4.874394
#> BV11-002-H_S58_MERGED.fastq.gz       354 13.982973  5.593423
#> BV11-002-M_S59_MERGED.fastq.gz        81 11.826913  7.395115
#> C21-NV1-B_S62_MERGED.fastq.gz         41 11.679756  5.974926
#> C21-NV1-H_S63_MERGED.fastq.gz        137  6.759908  2.621419
#> C21-NV1-M_S64_MERGED.fastq.gz         11  2.998873  1.806027
#> CB8-019-B_S69_MERGED.fastq.gz         25 16.491314 11.137183
#> CB8-019-H_S70_MERGED.fastq.gz         42  1.403082  1.109186
#> CB8-019-M_S71_MERGED.fastq.gz         24  4.000653  2.060547
#> D18-003-B_S78_MERGED.fastq.gz        136  9.725869  5.430373
#> D18-003-H_S79_MERGED.fastq.gz        111  2.672959  1.478804
#> D18-003-M_S80_MERGED.fastq.gz        127 10.072594  4.508733
#> D61-010-B_S82_MERGED.fastq.gz         32  1.527398  1.214526
#> D9-027-B_S83_MERGED.fastq.gz         112 19.748280  7.601955
#> D9-027-H_S84_MERGED.fastq.gz         118 12.838475  6.679902
#> D9-027-M_S85_MERGED.fastq.gz         151 24.272838  9.615561
#> DJ2-008-B_S87_MERGED.fastq.gz         49  8.487519  4.507450
#> DJ2-008-H_S88_MERGED.fastq.gz         44 18.454127 12.060999
#> DJ2-008-M_S89_MERGED.fastq.gz         77  5.650940  2.294316
#> DS1-ABM002-B_S91_MERGED.fastq.gz      91  7.468688  3.083528
#> DS1-ABM002-H_S92_MERGED.fastq.gz      46  3.468363  2.270055
#> DS1-ABM002-M_S93_MERGED.fastq.gz     156 18.696722 10.712748
#> DU3-045-B_S94_MERGED.fastq.gz        156 15.523102  7.153797
#> DY5-004-B_S96_MERGED.fastq.gz        100  7.340623  3.563105
#> DY5-004-H_S97_MERGED.fastq.gz         33  2.665689  1.664315
#> DY5-004-M_S98_MERGED.fastq.gz         20  4.771103  2.198900
#> E9-009-B_S100_MERGED.fastq.gz         31 24.852439 20.919368
#> E9-009-H_S101_MERGED.fastq.gz         44 20.687574 14.184711
#> E9-009-M_S102_MERGED.fastq.gz         32  2.781402  1.586266
#> EC2-013-B_S104_MERGED.fastq.gz       136  2.975489  1.535085
#> F7-015-M_S106_MERGED.fastq.gz         62  6.895465  3.522611
#> FOMES19-H_S108_MERGED.fastq.gz        33  1.576468  1.197767
#> FOMES19-M_S109_MERGED.fastq.gz        56  1.490947  1.160940
#> H10-018-M_S110_MERGED.fastq.gz         7  1.016148  1.003876
#> H24-NVABM1-H_S111_MERGED.fastq.gz     84 15.206040  7.015326
#> J18-004-B_S114_MERGED.fastq.gz        37 18.197369 11.411446
#> J18-004-H_S115_MERGED.fastq.gz        58  5.396765  3.287197
#> J18-004-M_S116_MERGED.fastq.gz        40 14.816708  9.053972
#> K18-002-H_S117_MERGED.fastq.gz       105  1.643799  1.171410
#> L19X-B_S119_MERGED.fastq.gz           70 20.115810 12.781959
#> L19X-H_S120_MERGED.fastq.gz          150 14.431518  7.461234
#> L19X-M_S121_MERGED.fastq.gz          127  9.432701  5.265973
#> L23-002-B_S122_MERGED.fastq.gz       181 18.608243  6.474391
#> L23-002-H_S123_MERGED.fastq.gz        62 10.618548  6.039820
#> L23-002-M_S124_MERGED.fastq.gz        33  3.282151  2.470113
#> N19X-B_S126_MERGED.fastq.gz           62 18.951406  8.222660
#> N19X-H_S127_MERGED.fastq.gz          257 19.214287  7.734104
#> N19X-M_S128_MERGED.fastq.gz          197 23.231017 11.207448
#> N22-001-B_S129_MERGED.fastq.gz         4  2.790604  2.314286
#> N23-002-B_S130_MERGED.fastq.gz       203  1.676570  1.165846
#> N23-002-H_S131_MERGED.fastq.gz        95 13.859756  7.945518
#> N23-002-M_S132_MERGED.fastq.gz         3  2.021246  1.750455
#> NVABM-0163-H_S135_MERGED.fastq.gz     91  9.966258  3.981573
#> NVABM0244-M_S137_MERGED.fastq.gz      58 12.056354  7.958510
#> O20-X-B_S139_MERGED.fastq.gz          31  4.989872  2.620331
#> O20-X-H_S140_MERGED.fastq.gz         184  5.393170  3.187427
#> O20-X-M_S141_MERGED.fastq.gz         146 15.558434  7.661789
#> O21-007-B_S142_MERGED.fastq.gz        22 15.238041  9.620712
#> O21-007-H_S143_MERGED.fastq.gz        20 16.559239 13.630903
#> O21-007-M_S144_MERGED.fastq.gz         5  2.742874  1.945632
#> O24-003-B_S145_MERGED.fastq.gz        38  4.863630  2.448777
#> O24-003-H_S146_MERGED.fastq.gz         9  7.985933  7.323964
#> O24-003-M_S147_MERGED.fastq.gz        30 14.554665  6.690734
#> O26-004-B_S148_MERGED.fastq.gz        61  3.237434  2.148826
#> O26-004-H_S149_MERGED.fastq.gz        25  1.243545  1.069097
#> O26-004-M_S150_MERGED.fastq.gz         9  2.455054  1.827384
#> O9-005-B_S152_MERGED.fastq.gz         62  4.870808  3.730590
#> P19-023-M_S153_MERGED.fastq.gz       149 11.756392  6.208680
#> P27-015-M_S154_MERGED.fastq.gz        73 12.110769  5.670984
#> Q27-ABM003-B_S156_MERGED.fastq.gz     64 16.985648  9.901538
#> R28-008-B_S158_MERGED.fastq.gz        53  5.060995  2.395426
#> R28-008-H_S159_MERGED.fastq.gz        24 16.126905 11.021638
#> R28-008-M_S160_MERGED.fastq.gz        31 24.605857 20.753838
#> T28-ABM602-B_S162_MERGED.fastq.gz    109  2.974004  1.727238
#> W26-001-B_S165_MERGED.fastq.gz        41  2.975811  1.773438
#> W26-001-H_S166_MERGED.fastq.gz        33 14.122368  6.129974
#> W26-001-M_S167_MERGED.fastq.gz        29 13.346043  6.752307
#> W9-025-M_S169_MERGED.fastq.gz         61  2.740403  1.529183
#> X24-009-B_S170_MERGED.fastq.gz        68  4.903310  2.590158
#> X24-009-H_S171_MERGED.fastq.gz        18  4.826186  4.353318
#> X24-009-M_S172_MERGED.fastq.gz        40  4.165627  3.655719
#> X29-004-B_S174_MERGED.fastq.gz        56  2.752176  1.605598
#> X29-004-H_S175_MERGED.fastq.gz       130 13.496186  6.609847
#> X29-004-M_S176_MERGED.fastq.gz       137  8.026196  4.372911
#> Y21-ABM484-H_S177_MERGED.fastq.gz     68  8.065746  4.665594
#> Y28-002-B_S178_MERGED.fastq.gz        46 16.108119  8.161933
#> Y28-002-H_S179_MERGED.fastq.gz        32 16.067341  8.609489
#> Y28-002-M_S180_MERGED.fastq.gz        18  4.809487  2.719143
#> Y29-007-B_S181_MERGED.fastq.gz        54 14.321086  6.100876
#> Y29-007-H_S182_MERGED.fastq.gz        26 13.893883  7.930828
#> Y29-007-M_S183_MERGED.fastq.gz        23 11.756051  7.119417
#> Y31-ABM484-B_S184_MERGED.fastq.gz    110  8.216186  5.130933
#> Z29-001-H_S185_MERGED.fastq.gz        85  7.813707  4.468417
#> Z30-ABM560-M_S187_MERGED.fastq.gz     80  5.659483  3.174734
#>                                                                   X
#> A10-005-B_S188_MERGED.fastq.gz       A10-005-B_S188_MERGED.fastq.gz
#> A10-005-H_S189_MERGED.fastq.gz       A10-005-H_S189_MERGED.fastq.gz
#> A10-005-M_S190_MERGED.fastq.gz       A10-005-M_S190_MERGED.fastq.gz
#> A12-007-B_S2_MERGED.fastq.gz           A12-007-B_S2_MERGED.fastq.gz
#> AB29-ABMX-H_S6_MERGED.fastq.gz       AB29-ABMX-H_S6_MERGED.fastq.gz
#> AD26-005-B_S9_MERGED.fastq.gz         AD26-005-B_S9_MERGED.fastq.gz
#> AD26-005-H_S10_MERGED.fastq.gz       AD26-005-H_S10_MERGED.fastq.gz
#> AD26-005-M_S11_MERGED.fastq.gz       AD26-005-M_S11_MERGED.fastq.gz
#> AD30-ABMX-M_S12_MERGED.fastq.gz     AD30-ABMX-M_S12_MERGED.fastq.gz
#> AD32-007-M_S13_MERGED.fastq.gz       AD32-007-M_S13_MERGED.fastq.gz
#> ADABM30X-B_S14_MERGED.fastq.gz       ADABM30X-B_S14_MERGED.fastq.gz
#> ADABM30X-H_S15_MERGED.fastq.gz       ADABM30X-H_S15_MERGED.fastq.gz
#> ADABM30X-M_S16_MERGED.fastq.gz       ADABM30X-M_S16_MERGED.fastq.gz
#> B18-006-B_S19_MERGED.fastq.gz         B18-006-B_S19_MERGED.fastq.gz
#> BA17-050-B_S21_MERGED.fastq.gz       BA17-050-B_S21_MERGED.fastq.gz
#> BB19-006-H_S22_MERGED.fastq.gz       BB19-006-H_S22_MERGED.fastq.gz
#> BB6-019-B_S23_MERGED.fastq.gz         BB6-019-B_S23_MERGED.fastq.gz
#> BB6-019-H_S24_MERGED.fastq.gz         BB6-019-H_S24_MERGED.fastq.gz
#> BB6-019-M_S25_MERGED.fastq.gz         BB6-019-M_S25_MERGED.fastq.gz
#> BE9-006-B_S27_MERGED.fastq.gz         BE9-006-B_S27_MERGED.fastq.gz
#> BE9-006-H_S28_MERGED.fastq.gz         BE9-006-H_S28_MERGED.fastq.gz
#> BE9-006-M_S29_MERGED.fastq.gz         BE9-006-M_S29_MERGED.fastq.gz
#> BG7-010-B_S30_MERGED.fastq.gz         BG7-010-B_S30_MERGED.fastq.gz
#> BG7-010-H_S31_MERGED.fastq.gz         BG7-010-H_S31_MERGED.fastq.gz
#> BG7-010-M_S32_MERGED.fastq.gz         BG7-010-M_S32_MERGED.fastq.gz
#> BJ17-007-M_S34_MERGED.fastq.gz       BJ17-007-M_S34_MERGED.fastq.gz
#> BL7-006-B_S36_MERGED.fastq.gz         BL7-006-B_S36_MERGED.fastq.gz
#> BL7-006-H_S37_MERGED.fastq.gz         BL7-006-H_S37_MERGED.fastq.gz
#> BL7-006-M_S38_MERGED.fastq.gz         BL7-006-M_S38_MERGED.fastq.gz
#> BP11-001-B_S43_MERGED.fastq.gz       BP11-001-B_S43_MERGED.fastq.gz
#> BP11-001-H_S44_MERGED.fastq.gz       BP11-001-H_S44_MERGED.fastq.gz
#> BP11-001-M_S45_MERGED.fastq.gz       BP11-001-M_S45_MERGED.fastq.gz
#> BP12-025-B_S46_MERGED.fastq.gz       BP12-025-B_S46_MERGED.fastq.gz
#> BQ4-018-B_S49_MERGED.fastq.gz         BQ4-018-B_S49_MERGED.fastq.gz
#> BQ4-018-H_S50_MERGED.fastq.gz         BQ4-018-H_S50_MERGED.fastq.gz
#> BQ4-018-M_S51_MERGED.fastq.gz         BQ4-018-M_S51_MERGED.fastq.gz
#> BT-006-M_S55_MERGED.fastq.gz           BT-006-M_S55_MERGED.fastq.gz
#> BV11-002-B_S57_MERGED.fastq.gz       BV11-002-B_S57_MERGED.fastq.gz
#> BV11-002-H_S58_MERGED.fastq.gz       BV11-002-H_S58_MERGED.fastq.gz
#> BV11-002-M_S59_MERGED.fastq.gz       BV11-002-M_S59_MERGED.fastq.gz
#> C21-NV1-B_S62_MERGED.fastq.gz         C21-NV1-B_S62_MERGED.fastq.gz
#> C21-NV1-H_S63_MERGED.fastq.gz         C21-NV1-H_S63_MERGED.fastq.gz
#> C21-NV1-M_S64_MERGED.fastq.gz         C21-NV1-M_S64_MERGED.fastq.gz
#> CB8-019-B_S69_MERGED.fastq.gz         CB8-019-B_S69_MERGED.fastq.gz
#> CB8-019-H_S70_MERGED.fastq.gz         CB8-019-H_S70_MERGED.fastq.gz
#> CB8-019-M_S71_MERGED.fastq.gz         CB8-019-M_S71_MERGED.fastq.gz
#> D18-003-B_S78_MERGED.fastq.gz         D18-003-B_S78_MERGED.fastq.gz
#> D18-003-H_S79_MERGED.fastq.gz         D18-003-H_S79_MERGED.fastq.gz
#> D18-003-M_S80_MERGED.fastq.gz         D18-003-M_S80_MERGED.fastq.gz
#> D61-010-B_S82_MERGED.fastq.gz         D61-010-B_S82_MERGED.fastq.gz
#> D9-027-B_S83_MERGED.fastq.gz           D9-027-B_S83_MERGED.fastq.gz
#> D9-027-H_S84_MERGED.fastq.gz           D9-027-H_S84_MERGED.fastq.gz
#> D9-027-M_S85_MERGED.fastq.gz           D9-027-M_S85_MERGED.fastq.gz
#> DJ2-008-B_S87_MERGED.fastq.gz         DJ2-008-B_S87_MERGED.fastq.gz
#> DJ2-008-H_S88_MERGED.fastq.gz         DJ2-008-H_S88_MERGED.fastq.gz
#> DJ2-008-M_S89_MERGED.fastq.gz         DJ2-008-M_S89_MERGED.fastq.gz
#> DS1-ABM002-B_S91_MERGED.fastq.gz   DS1-ABM002-B_S91_MERGED.fastq.gz
#> DS1-ABM002-H_S92_MERGED.fastq.gz   DS1-ABM002-H_S92_MERGED.fastq.gz
#> DS1-ABM002-M_S93_MERGED.fastq.gz   DS1-ABM002-M_S93_MERGED.fastq.gz
#> DU3-045-B_S94_MERGED.fastq.gz         DU3-045-B_S94_MERGED.fastq.gz
#> DY5-004-B_S96_MERGED.fastq.gz         DY5-004-B_S96_MERGED.fastq.gz
#> DY5-004-H_S97_MERGED.fastq.gz         DY5-004-H_S97_MERGED.fastq.gz
#> DY5-004-M_S98_MERGED.fastq.gz         DY5-004-M_S98_MERGED.fastq.gz
#> E9-009-B_S100_MERGED.fastq.gz         E9-009-B_S100_MERGED.fastq.gz
#> E9-009-H_S101_MERGED.fastq.gz         E9-009-H_S101_MERGED.fastq.gz
#> E9-009-M_S102_MERGED.fastq.gz         E9-009-M_S102_MERGED.fastq.gz
#> EC2-013-B_S104_MERGED.fastq.gz       EC2-013-B_S104_MERGED.fastq.gz
#> F7-015-M_S106_MERGED.fastq.gz         F7-015-M_S106_MERGED.fastq.gz
#> FOMES19-H_S108_MERGED.fastq.gz       FOMES19-H_S108_MERGED.fastq.gz
#> FOMES19-M_S109_MERGED.fastq.gz       FOMES19-M_S109_MERGED.fastq.gz
#> H10-018-M_S110_MERGED.fastq.gz       H10-018-M_S110_MERGED.fastq.gz
#> H24-NVABM1-H_S111_MERGED.fastq.gz H24-NVABM1-H_S111_MERGED.fastq.gz
#> J18-004-B_S114_MERGED.fastq.gz       J18-004-B_S114_MERGED.fastq.gz
#> J18-004-H_S115_MERGED.fastq.gz       J18-004-H_S115_MERGED.fastq.gz
#> J18-004-M_S116_MERGED.fastq.gz       J18-004-M_S116_MERGED.fastq.gz
#> K18-002-H_S117_MERGED.fastq.gz       K18-002-H_S117_MERGED.fastq.gz
#> L19X-B_S119_MERGED.fastq.gz             L19X-B_S119_MERGED.fastq.gz
#> L19X-H_S120_MERGED.fastq.gz             L19X-H_S120_MERGED.fastq.gz
#> L19X-M_S121_MERGED.fastq.gz             L19X-M_S121_MERGED.fastq.gz
#> L23-002-B_S122_MERGED.fastq.gz       L23-002-B_S122_MERGED.fastq.gz
#> L23-002-H_S123_MERGED.fastq.gz       L23-002-H_S123_MERGED.fastq.gz
#> L23-002-M_S124_MERGED.fastq.gz       L23-002-M_S124_MERGED.fastq.gz
#> N19X-B_S126_MERGED.fastq.gz             N19X-B_S126_MERGED.fastq.gz
#> N19X-H_S127_MERGED.fastq.gz             N19X-H_S127_MERGED.fastq.gz
#> N19X-M_S128_MERGED.fastq.gz             N19X-M_S128_MERGED.fastq.gz
#> N22-001-B_S129_MERGED.fastq.gz       N22-001-B_S129_MERGED.fastq.gz
#> N23-002-B_S130_MERGED.fastq.gz       N23-002-B_S130_MERGED.fastq.gz
#> N23-002-H_S131_MERGED.fastq.gz       N23-002-H_S131_MERGED.fastq.gz
#> N23-002-M_S132_MERGED.fastq.gz       N23-002-M_S132_MERGED.fastq.gz
#> NVABM-0163-H_S135_MERGED.fastq.gz NVABM-0163-H_S135_MERGED.fastq.gz
#> NVABM0244-M_S137_MERGED.fastq.gz   NVABM0244-M_S137_MERGED.fastq.gz
#> O20-X-B_S139_MERGED.fastq.gz           O20-X-B_S139_MERGED.fastq.gz
#> O20-X-H_S140_MERGED.fastq.gz           O20-X-H_S140_MERGED.fastq.gz
#> O20-X-M_S141_MERGED.fastq.gz           O20-X-M_S141_MERGED.fastq.gz
#> O21-007-B_S142_MERGED.fastq.gz       O21-007-B_S142_MERGED.fastq.gz
#> O21-007-H_S143_MERGED.fastq.gz       O21-007-H_S143_MERGED.fastq.gz
#> O21-007-M_S144_MERGED.fastq.gz       O21-007-M_S144_MERGED.fastq.gz
#> O24-003-B_S145_MERGED.fastq.gz       O24-003-B_S145_MERGED.fastq.gz
#> O24-003-H_S146_MERGED.fastq.gz       O24-003-H_S146_MERGED.fastq.gz
#> O24-003-M_S147_MERGED.fastq.gz       O24-003-M_S147_MERGED.fastq.gz
#> O26-004-B_S148_MERGED.fastq.gz       O26-004-B_S148_MERGED.fastq.gz
#> O26-004-H_S149_MERGED.fastq.gz       O26-004-H_S149_MERGED.fastq.gz
#> O26-004-M_S150_MERGED.fastq.gz       O26-004-M_S150_MERGED.fastq.gz
#> O9-005-B_S152_MERGED.fastq.gz         O9-005-B_S152_MERGED.fastq.gz
#> P19-023-M_S153_MERGED.fastq.gz       P19-023-M_S153_MERGED.fastq.gz
#> P27-015-M_S154_MERGED.fastq.gz       P27-015-M_S154_MERGED.fastq.gz
#> Q27-ABM003-B_S156_MERGED.fastq.gz Q27-ABM003-B_S156_MERGED.fastq.gz
#> R28-008-B_S158_MERGED.fastq.gz       R28-008-B_S158_MERGED.fastq.gz
#> R28-008-H_S159_MERGED.fastq.gz       R28-008-H_S159_MERGED.fastq.gz
#> R28-008-M_S160_MERGED.fastq.gz       R28-008-M_S160_MERGED.fastq.gz
#> T28-ABM602-B_S162_MERGED.fastq.gz T28-ABM602-B_S162_MERGED.fastq.gz
#> W26-001-B_S165_MERGED.fastq.gz       W26-001-B_S165_MERGED.fastq.gz
#> W26-001-H_S166_MERGED.fastq.gz       W26-001-H_S166_MERGED.fastq.gz
#> W26-001-M_S167_MERGED.fastq.gz       W26-001-M_S167_MERGED.fastq.gz
#> W9-025-M_S169_MERGED.fastq.gz         W9-025-M_S169_MERGED.fastq.gz
#> X24-009-B_S170_MERGED.fastq.gz       X24-009-B_S170_MERGED.fastq.gz
#> X24-009-H_S171_MERGED.fastq.gz       X24-009-H_S171_MERGED.fastq.gz
#> X24-009-M_S172_MERGED.fastq.gz       X24-009-M_S172_MERGED.fastq.gz
#> X29-004-B_S174_MERGED.fastq.gz       X29-004-B_S174_MERGED.fastq.gz
#> X29-004-H_S175_MERGED.fastq.gz       X29-004-H_S175_MERGED.fastq.gz
#> X29-004-M_S176_MERGED.fastq.gz       X29-004-M_S176_MERGED.fastq.gz
#> Y21-ABM484-H_S177_MERGED.fastq.gz Y21-ABM484-H_S177_MERGED.fastq.gz
#> Y28-002-B_S178_MERGED.fastq.gz       Y28-002-B_S178_MERGED.fastq.gz
#> Y28-002-H_S179_MERGED.fastq.gz       Y28-002-H_S179_MERGED.fastq.gz
#> Y28-002-M_S180_MERGED.fastq.gz       Y28-002-M_S180_MERGED.fastq.gz
#> Y29-007-B_S181_MERGED.fastq.gz       Y29-007-B_S181_MERGED.fastq.gz
#> Y29-007-H_S182_MERGED.fastq.gz       Y29-007-H_S182_MERGED.fastq.gz
#> Y29-007-M_S183_MERGED.fastq.gz       Y29-007-M_S183_MERGED.fastq.gz
#> Y31-ABM484-B_S184_MERGED.fastq.gz Y31-ABM484-B_S184_MERGED.fastq.gz
#> Z29-001-H_S185_MERGED.fastq.gz       Z29-001-H_S185_MERGED.fastq.gz
#> Z30-ABM560-M_S187_MERGED.fastq.gz Z30-ABM560-M_S187_MERGED.fastq.gz
#>                                        Sample_names   Tree_name Sample_id
#> A10-005-B_S188_MERGED.fastq.gz       A10-005-B_S188     A10-005       188
#> A10-005-H_S189_MERGED.fastq.gz       A10-005-H_S189     A10-005       189
#> A10-005-M_S190_MERGED.fastq.gz       A10-005-M_S190     A10-005       190
#> A12-007-B_S2_MERGED.fastq.gz           A12-007-B_S2     A12-007         2
#> AB29-ABMX-H_S6_MERGED.fastq.gz       AB29-ABMX-H_S6  AB29-abm-X         6
#> AD26-005-B_S9_MERGED.fastq.gz         AD26-005-B_S9    AD26-005         9
#> AD26-005-H_S10_MERGED.fastq.gz       AD26-005-H_S10    AD26-005        10
#> AD26-005-M_S11_MERGED.fastq.gz       AD26-005-M_S11    AD26-005        11
#> AD30-ABMX-M_S12_MERGED.fastq.gz     AD30-ABMX-M_S12  AD30-abm-X        12
#> AD32-007-M_S13_MERGED.fastq.gz       AD32-007-M_S13    AD32-007        13
#> ADABM30X-B_S14_MERGED.fastq.gz       ADABM30X-B_S14  AD30-abm-X        14
#> ADABM30X-H_S15_MERGED.fastq.gz       ADABM30X-H_S15  AD30-abm-X        15
#> ADABM30X-M_S16_MERGED.fastq.gz       ADABM30X-M_S16  AD30-abm-X        16
#> B18-006-B_S19_MERGED.fastq.gz         B18-006-B_S19     B18-006        19
#> BA17-050-B_S21_MERGED.fastq.gz       BA17-050-B_S21    BA17-050        21
#> BB19-006-H_S22_MERGED.fastq.gz       BB19-006-H_S22    BB19-006        22
#> BB6-019-B_S23_MERGED.fastq.gz         BB6-019-B_S23     BB6-019        23
#> BB6-019-H_S24_MERGED.fastq.gz         BB6-019-H_S24     BB6-019        24
#> BB6-019-M_S25_MERGED.fastq.gz         BB6-019-M_S25     BB6-019        25
#> BE9-006-B_S27_MERGED.fastq.gz         BE9-006-B_S27     BE9-006        27
#> BE9-006-H_S28_MERGED.fastq.gz         BE9-006-H_S28     BE9-006        28
#> BE9-006-M_S29_MERGED.fastq.gz         BE9-006-M_S29     BE9-006        29
#> BG7-010-B_S30_MERGED.fastq.gz         BG7-010-B_S30     BG7-010        30
#> BG7-010-H_S31_MERGED.fastq.gz         BG7-010-H_S31     BG7-010        31
#> BG7-010-M_S32_MERGED.fastq.gz         BG7-010-M_S32     BG7-010        32
#> BJ17-007-M_S34_MERGED.fastq.gz       BJ17-007-M_S34    BJ17-007        34
#> BL7-006-B_S36_MERGED.fastq.gz         BL7-006-B_S36     BL7-006        36
#> BL7-006-H_S37_MERGED.fastq.gz         BL7-006-H_S37     BL7-006        37
#> BL7-006-M_S38_MERGED.fastq.gz         BL7-006-M_S38     BL7-006        38
#> BP11-001-B_S43_MERGED.fastq.gz       BP11-001-B_S43    BP11-001        43
#> BP11-001-H_S44_MERGED.fastq.gz       BP11-001-H_S44    BP11-001        44
#> BP11-001-M_S45_MERGED.fastq.gz       BP11-001-M_S45    BP11-001        45
#> BP12-025-B_S46_MERGED.fastq.gz       BP12-025-B_S46    BP12-025        46
#> BQ4-018-B_S49_MERGED.fastq.gz         BQ4-018-B_S49     BQ4-018        49
#> BQ4-018-H_S50_MERGED.fastq.gz         BQ4-018-H_S50     BQ4-018        50
#> BQ4-018-M_S51_MERGED.fastq.gz         BQ4-018-M_S51     BQ4-018        51
#> BT-006-M_S55_MERGED.fastq.gz           BT-006-M_S55     BT7-006        55
#> BV11-002-B_S57_MERGED.fastq.gz       BV11-002-B_S57    BU11-002        57
#> BV11-002-H_S58_MERGED.fastq.gz       BV11-002-H_S58    BU11-002        58
#> BV11-002-M_S59_MERGED.fastq.gz       BV11-002-M_S59    BU11-002        59
#> C21-NV1-B_S62_MERGED.fastq.gz         C21-NV1-B_S62    C21-nv-1        62
#> C21-NV1-H_S63_MERGED.fastq.gz         C21-NV1-H_S63    C21-nv-1        63
#> C21-NV1-M_S64_MERGED.fastq.gz         C21-NV1-M_S64    C21-nv-1        64
#> CB8-019-B_S69_MERGED.fastq.gz         CB8-019-B_S69     CB8-019        69
#> CB8-019-H_S70_MERGED.fastq.gz         CB8-019-H_S70     CB8-019        70
#> CB8-019-M_S71_MERGED.fastq.gz         CB8-019-M_S71     CB8-019        71
#> D18-003-B_S78_MERGED.fastq.gz         D18-003-B_S78     D18-003        78
#> D18-003-H_S79_MERGED.fastq.gz         D18-003-H_S79     D18-003        79
#> D18-003-M_S80_MERGED.fastq.gz         D18-003-M_S80     D18-003        80
#> D61-010-B_S82_MERGED.fastq.gz         D61-010-B_S82     DG1-010        82
#> D9-027-B_S83_MERGED.fastq.gz           D9-027-B_S83      D9-027        83
#> D9-027-H_S84_MERGED.fastq.gz           D9-027-H_S84      D9-027        84
#> D9-027-M_S85_MERGED.fastq.gz           D9-027-M_S85      D9-027        85
#> DJ2-008-B_S87_MERGED.fastq.gz         DJ2-008-B_S87     DJ2-008        87
#> DJ2-008-H_S88_MERGED.fastq.gz         DJ2-008-H_S88     DJ2-008        88
#> DJ2-008-M_S89_MERGED.fastq.gz         DJ2-008-M_S89     DJ2-008        89
#> DS1-ABM002-B_S91_MERGED.fastq.gz   DS1-ABM002-B_S91   DS1abm002        91
#> DS1-ABM002-H_S92_MERGED.fastq.gz   DS1-ABM002-H_S92   DS1abm002        92
#> DS1-ABM002-M_S93_MERGED.fastq.gz   DS1-ABM002-M_S93   DS1abm002        93
#> DU3-045-B_S94_MERGED.fastq.gz         DU3-045-B_S94     DU3-045        94
#> DY5-004-B_S96_MERGED.fastq.gz         DY5-004-B_S96     DY5-004        96
#> DY5-004-H_S97_MERGED.fastq.gz         DY5-004-H_S97     DY5-004        97
#> DY5-004-M_S98_MERGED.fastq.gz         DY5-004-M_S98     DY5-004        98
#> E9-009-B_S100_MERGED.fastq.gz         E9-009-B_S100      E9-009       100
#> E9-009-H_S101_MERGED.fastq.gz         E9-009-H_S101      E9-009       101
#> E9-009-M_S102_MERGED.fastq.gz         E9-009-M_S102      E9-009       102
#> EC2-013-B_S104_MERGED.fastq.gz       EC2-013-B_S104     EC2-013       104
#> F7-015-M_S106_MERGED.fastq.gz         F7-015-M_S106      F7-015       106
#> FOMES19-H_S108_MERGED.fastq.gz       FOMES19-H_S108     FOMES19       108
#> FOMES19-M_S109_MERGED.fastq.gz       FOMES19-M_S109     FOMES19       109
#> H10-018-M_S110_MERGED.fastq.gz       H10-018-M_S110     H10-018       110
#> H24-NVABM1-H_S111_MERGED.fastq.gz H24-NVABM1-H_S111 H24-nvabm-1       111
#> J18-004-B_S114_MERGED.fastq.gz       J18-004-B_S114     J18-004       114
#> J18-004-H_S115_MERGED.fastq.gz       J18-004-H_S115     J18-004       115
#> J18-004-M_S116_MERGED.fastq.gz       J18-004-M_S116     J18-004       116
#> K18-002-H_S117_MERGED.fastq.gz       K18-002-H_S117     K18-002       117
#> L19X-B_S119_MERGED.fastq.gz             L19X-B_S119       L19-X       119
#> L19X-H_S120_MERGED.fastq.gz             L19X-H_S120       L19-X       120
#> L19X-M_S121_MERGED.fastq.gz             L19X-M_S121       L19-X       121
#> L23-002-B_S122_MERGED.fastq.gz       L23-002-B_S122     L23-002       122
#> L23-002-H_S123_MERGED.fastq.gz       L23-002-H_S123     L23-002       123
#> L23-002-M_S124_MERGED.fastq.gz       L23-002-M_S124     L23-002       124
#> N19X-B_S126_MERGED.fastq.gz             N19X-B_S126       N19-X       126
#> N19X-H_S127_MERGED.fastq.gz             N19X-H_S127       N19-X       127
#> N19X-M_S128_MERGED.fastq.gz             N19X-M_S128       N19-X       128
#> N22-001-B_S129_MERGED.fastq.gz       N22-001-B_S129     M22-001       129
#> N23-002-B_S130_MERGED.fastq.gz       N23-002-B_S130     N23-002       130
#> N23-002-H_S131_MERGED.fastq.gz       N23-002-H_S131     N23-002       131
#> N23-002-M_S132_MERGED.fastq.gz       N23-002-M_S132     N23-002       132
#> NVABM-0163-H_S135_MERGED.fastq.gz NVABM-0163-H_S135   nvabm0163       135
#> NVABM0244-M_S137_MERGED.fastq.gz   NVABM0244-M_S137   nvabm0244       137
#> O20-X-B_S139_MERGED.fastq.gz           O20-X-B_S139       O20-X       139
#> O20-X-H_S140_MERGED.fastq.gz           O20-X-H_S140       O20-X       140
#> O20-X-M_S141_MERGED.fastq.gz           O20-X-M_S141       O20-X       141
#> O21-007-B_S142_MERGED.fastq.gz       O21-007-B_S142     O21-007       142
#> O21-007-H_S143_MERGED.fastq.gz       O21-007-H_S143     O21-007       143
#> O21-007-M_S144_MERGED.fastq.gz       O21-007-M_S144     O21-007       144
#> O24-003-B_S145_MERGED.fastq.gz       O24-003-B_S145     O24-003       145
#> O24-003-H_S146_MERGED.fastq.gz       O24-003-H_S146     O24-003       146
#> O24-003-M_S147_MERGED.fastq.gz       O24-003-M_S147     O24-003       147
#> O26-004-B_S148_MERGED.fastq.gz       O26-004-B_S148     O26-004       148
#> O26-004-H_S149_MERGED.fastq.gz       O26-004-H_S149     O26-004       149
#> O26-004-M_S150_MERGED.fastq.gz       O26-004-M_S150     O26-004       150
#> O9-005-B_S152_MERGED.fastq.gz         O9-005-B_S152      C9-005       152
#> P19-023-M_S153_MERGED.fastq.gz       P19-023-M_S153     P19-023       153
#> P27-015-M_S154_MERGED.fastq.gz       P27-015-M_S154     P27-015       154
#> Q27-ABM003-B_S156_MERGED.fastq.gz Q27-ABM003-B_S156  Q27-abm003       156
#> R28-008-B_S158_MERGED.fastq.gz       R28-008-B_S158     R28-008       158
#> R28-008-H_S159_MERGED.fastq.gz       R28-008-H_S159     R28-008       159
#> R28-008-M_S160_MERGED.fastq.gz       R28-008-M_S160     R28-008       160
#> T28-ABM602-B_S162_MERGED.fastq.gz T28-ABM602-B_S162  T28-abm602       162
#> W26-001-B_S165_MERGED.fastq.gz       W26-001-B_S165     W26-001       165
#> W26-001-H_S166_MERGED.fastq.gz       W26-001-H_S166     W26-001       166
#> W26-001-M_S167_MERGED.fastq.gz       W26-001-M_S167     W26-001       167
#> W9-025-M_S169_MERGED.fastq.gz         W9-025-M_S169     W29-025       169
#> X24-009-B_S170_MERGED.fastq.gz       X24-009-B_S170     X24-009       170
#> X24-009-H_S171_MERGED.fastq.gz       X24-009-H_S171     X24-009       171
#> X24-009-M_S172_MERGED.fastq.gz       X24-009-M_S172     X24-009       172
#> X29-004-B_S174_MERGED.fastq.gz       X29-004-B_S174     X29-004       174
#> X29-004-H_S175_MERGED.fastq.gz       X29-004-H_S175     X29-004       175
#> X29-004-M_S176_MERGED.fastq.gz       X29-004-M_S176     X29-004       176
#> Y21-ABM484-H_S177_MERGED.fastq.gz Y21-ABM484-H_S177  Y31-abm484       177
#> Y28-002-B_S178_MERGED.fastq.gz       Y28-002-B_S178     Y28-002       178
#> Y28-002-H_S179_MERGED.fastq.gz       Y28-002-H_S179     Y28-002       179
#> Y28-002-M_S180_MERGED.fastq.gz       Y28-002-M_S180     Y28-002       180
#> Y29-007-B_S181_MERGED.fastq.gz       Y29-007-B_S181     Y29-007       181
#> Y29-007-H_S182_MERGED.fastq.gz       Y29-007-H_S182     Y29-007       182
#> Y29-007-M_S183_MERGED.fastq.gz       Y29-007-M_S183     Y29-007       183
#> Y31-ABM484-B_S184_MERGED.fastq.gz Y31-ABM484-B_S184  Y31-abm484       184
#> Z29-001-H_S185_MERGED.fastq.gz       Z29-001-H_S185     Z29-001       185
#> Z30-ABM560-M_S187_MERGED.fastq.gz Z30-ABM560-M_S187  Z30-abm560       187
#>                                   Height Diameter Time
#> A10-005-B_S188_MERGED.fastq.gz       Low       52   15
#> A10-005-H_S189_MERGED.fastq.gz      High       52   15
#> A10-005-M_S190_MERGED.fastq.gz    Middle       52   15
#> A12-007-B_S2_MERGED.fastq.gz         Low     28,4    0
#> AB29-ABMX-H_S6_MERGED.fastq.gz      High       99    5
#> AD26-005-B_S9_MERGED.fastq.gz        Low    115,5   15
#> AD26-005-H_S10_MERGED.fastq.gz      High    115,5   15
#> AD26-005-M_S11_MERGED.fastq.gz    Middle    115,5   15
#> AD30-ABMX-M_S12_MERGED.fastq.gz   Middle        -    5
#> AD32-007-M_S13_MERGED.fastq.gz    Middle       52   NA
#> ADABM30X-B_S14_MERGED.fastq.gz       Low        -    5
#> ADABM30X-H_S15_MERGED.fastq.gz      High        -    5
#> ADABM30X-M_S16_MERGED.fastq.gz    Middle        -    5
#> B18-006-B_S19_MERGED.fastq.gz        Low     32,5   10
#> BA17-050-B_S21_MERGED.fastq.gz       Low       30   10
#> BB19-006-H_S22_MERGED.fastq.gz      High     37,5   15
#> BB6-019-B_S23_MERGED.fastq.gz        Low     33,8    5
#> BB6-019-H_S24_MERGED.fastq.gz       High     33,8    5
#> BB6-019-M_S25_MERGED.fastq.gz     Middle     33,8    5
#> BE9-006-B_S27_MERGED.fastq.gz        Low     23,1    0
#> BE9-006-H_S28_MERGED.fastq.gz       High     23,1    0
#> BE9-006-M_S29_MERGED.fastq.gz     Middle     23,1    0
#> BG7-010-B_S30_MERGED.fastq.gz        Low       32    0
#> BG7-010-H_S31_MERGED.fastq.gz       High       32    0
#> BG7-010-M_S32_MERGED.fastq.gz     Middle       32    0
#> BJ17-007-M_S34_MERGED.fastq.gz    Middle       32    5
#> BL7-006-B_S36_MERGED.fastq.gz        Low     35,6    0
#> BL7-006-H_S37_MERGED.fastq.gz       High     35,6    0
#> BL7-006-M_S38_MERGED.fastq.gz     Middle     35,6    0
#> BP11-001-B_S43_MERGED.fastq.gz       Low     47,5   15
#> BP11-001-H_S44_MERGED.fastq.gz      High     47,5   15
#> BP11-001-M_S45_MERGED.fastq.gz    Middle     47,5   15
#> BP12-025-B_S46_MERGED.fastq.gz       Low     64,5   10
#> BQ4-018-B_S49_MERGED.fastq.gz        Low       34    0
#> BQ4-018-H_S50_MERGED.fastq.gz       High       34    0
#> BQ4-018-M_S51_MERGED.fastq.gz     Middle       34    0
#> BT-006-M_S55_MERGED.fastq.gz      Middle       41    0
#> BV11-002-B_S57_MERGED.fastq.gz       Low       33    5
#> BV11-002-H_S58_MERGED.fastq.gz      High       33    5
#> BV11-002-M_S59_MERGED.fastq.gz    Middle       33    5
#> C21-NV1-B_S62_MERGED.fastq.gz        Low       30    0
#> C21-NV1-H_S63_MERGED.fastq.gz       High       30    0
#> C21-NV1-M_S64_MERGED.fastq.gz     Middle       30    0
#> CB8-019-B_S69_MERGED.fastq.gz        Low     33,3    0
#> CB8-019-H_S70_MERGED.fastq.gz       High     33,3    0
#> CB8-019-M_S71_MERGED.fastq.gz     Middle     33,3    0
#> D18-003-B_S78_MERGED.fastq.gz        Low     45,5    5
#> D18-003-H_S79_MERGED.fastq.gz       High     45,5    5
#> D18-003-M_S80_MERGED.fastq.gz     Middle     45,5    5
#> D61-010-B_S82_MERGED.fastq.gz        Low     36,3   15
#> D9-027-B_S83_MERGED.fastq.gz         Low     18,9   NA
#> D9-027-H_S84_MERGED.fastq.gz        High     18,9   NA
#> D9-027-M_S85_MERGED.fastq.gz      Middle     18,9   NA
#> DJ2-008-B_S87_MERGED.fastq.gz        Low     76,4    0
#> DJ2-008-H_S88_MERGED.fastq.gz       High     76,4    0
#> DJ2-008-M_S89_MERGED.fastq.gz     Middle     76,4    0
#> DS1-ABM002-B_S91_MERGED.fastq.gz     Low       25   15
#> DS1-ABM002-H_S92_MERGED.fastq.gz    High       25   15
#> DS1-ABM002-M_S93_MERGED.fastq.gz  Middle       25   15
#> DU3-045-B_S94_MERGED.fastq.gz        Low       38   15
#> DY5-004-B_S96_MERGED.fastq.gz        Low     60,4    0
#> DY5-004-H_S97_MERGED.fastq.gz       High     60,4    0
#> DY5-004-M_S98_MERGED.fastq.gz     Middle     60,4    0
#> E9-009-B_S100_MERGED.fastq.gz        Low     29,5    0
#> E9-009-H_S101_MERGED.fastq.gz       High     29,5    0
#> E9-009-M_S102_MERGED.fastq.gz     Middle     29,5    0
#> EC2-013-B_S104_MERGED.fastq.gz       Low       52   15
#> F7-015-M_S106_MERGED.fastq.gz     Middle       50   10
#> FOMES19-H_S108_MERGED.fastq.gz      High     <NA>   NA
#> FOMES19-M_S109_MERGED.fastq.gz    Middle     <NA>   NA
#> H10-018-M_S110_MERGED.fastq.gz    Middle     88,5   15
#> H24-NVABM1-H_S111_MERGED.fastq.gz   High       30    5
#> J18-004-B_S114_MERGED.fastq.gz       Low     27,6    0
#> J18-004-H_S115_MERGED.fastq.gz      High     27,6    0
#> J18-004-M_S116_MERGED.fastq.gz    Middle     27,6    0
#> K18-002-H_S117_MERGED.fastq.gz      High     33,5   10
#> L19X-B_S119_MERGED.fastq.gz          Low       30   NA
#> L19X-H_S120_MERGED.fastq.gz         High       30   NA
#> L19X-M_S121_MERGED.fastq.gz       Middle       30   NA
#> L23-002-B_S122_MERGED.fastq.gz       Low     40,5   10
#> L23-002-H_S123_MERGED.fastq.gz      High     40,5   10
#> L23-002-M_S124_MERGED.fastq.gz    Middle     40,5   10
#> N19X-B_S126_MERGED.fastq.gz          Low     33,5   NA
#> N19X-H_S127_MERGED.fastq.gz         High     33,5   NA
#> N19X-M_S128_MERGED.fastq.gz       Middle     33,5   NA
#> N22-001-B_S129_MERGED.fastq.gz       Low       34    0
#> N23-002-B_S130_MERGED.fastq.gz       Low     21,5   15
#> N23-002-H_S131_MERGED.fastq.gz      High     21,5   15
#> N23-002-M_S132_MERGED.fastq.gz    Middle     21,5   15
#> NVABM-0163-H_S135_MERGED.fastq.gz   High       40   10
#> NVABM0244-M_S137_MERGED.fastq.gz  Middle       50   10
#> O20-X-B_S139_MERGED.fastq.gz         Low     <NA>   NA
#> O20-X-H_S140_MERGED.fastq.gz        High     <NA>   NA
#> O20-X-M_S141_MERGED.fastq.gz      Middle     <NA>   NA
#> O21-007-B_S142_MERGED.fastq.gz       Low       25    0
#> O21-007-H_S143_MERGED.fastq.gz      High       25    0
#> O21-007-M_S144_MERGED.fastq.gz    Middle       25    0
#> O24-003-B_S145_MERGED.fastq.gz       Low       75    0
#> O24-003-H_S146_MERGED.fastq.gz      High       75    0
#> O24-003-M_S147_MERGED.fastq.gz    Middle       75    0
#> O26-004-B_S148_MERGED.fastq.gz       Low     30,4    0
#> O26-004-H_S149_MERGED.fastq.gz      High     30,4    0
#> O26-004-M_S150_MERGED.fastq.gz    Middle     30,4    0
#> O9-005-B_S152_MERGED.fastq.gz        Low       39    5
#> P19-023-M_S153_MERGED.fastq.gz    Middle     18,4   NA
#> P27-015-M_S154_MERGED.fastq.gz    Middle       96    5
#> Q27-ABM003-B_S156_MERGED.fastq.gz    Low       10    5
#> R28-008-B_S158_MERGED.fastq.gz       Low     37,2    0
#> R28-008-H_S159_MERGED.fastq.gz      High     37,2    0
#> R28-008-M_S160_MERGED.fastq.gz    Middle     37,2    0
#> T28-ABM602-B_S162_MERGED.fastq.gz    Low       64    5
#> W26-001-B_S165_MERGED.fastq.gz       Low     71,7    0
#> W26-001-H_S166_MERGED.fastq.gz      High     71,7    0
#> W26-001-M_S167_MERGED.fastq.gz    Middle     71,7    0
#> W9-025-M_S169_MERGED.fastq.gz     Middle     71,5   NA
#> X24-009-B_S170_MERGED.fastq.gz       Low     38,8    0
#> X24-009-H_S171_MERGED.fastq.gz      High     38,8    0
#> X24-009-M_S172_MERGED.fastq.gz    Middle     38,8    0
#> X29-004-B_S174_MERGED.fastq.gz       Low     52,7   10
#> X29-004-H_S175_MERGED.fastq.gz      High     52,7   10
#> X29-004-M_S176_MERGED.fastq.gz    Middle     52,7   10
#> Y21-ABM484-H_S177_MERGED.fastq.gz   High       20    5
#> Y28-002-B_S178_MERGED.fastq.gz       Low     63,6   NA
#> Y28-002-H_S179_MERGED.fastq.gz      High     63,6   NA
#> Y28-002-M_S180_MERGED.fastq.gz    Middle     63,6   NA
#> Y29-007-B_S181_MERGED.fastq.gz       Low     68,6    0
#> Y29-007-H_S182_MERGED.fastq.gz      High     68,6    0
#> Y29-007-M_S183_MERGED.fastq.gz    Middle     68,6    0
#> Y31-ABM484-B_S184_MERGED.fastq.gz    Low       20    5
#> Z29-001-H_S185_MERGED.fastq.gz      High     35,2    5
#> Z30-ABM560-M_S187_MERGED.fastq.gz Middle       11    5
#>                                                          .sample_id
#> A10-005-B_S188_MERGED.fastq.gz       A10-005-B_S188_MERGED.fastq.gz
#> A10-005-H_S189_MERGED.fastq.gz       A10-005-H_S189_MERGED.fastq.gz
#> A10-005-M_S190_MERGED.fastq.gz       A10-005-M_S190_MERGED.fastq.gz
#> A12-007-B_S2_MERGED.fastq.gz           A12-007-B_S2_MERGED.fastq.gz
#> AB29-ABMX-H_S6_MERGED.fastq.gz       AB29-ABMX-H_S6_MERGED.fastq.gz
#> AD26-005-B_S9_MERGED.fastq.gz         AD26-005-B_S9_MERGED.fastq.gz
#> AD26-005-H_S10_MERGED.fastq.gz       AD26-005-H_S10_MERGED.fastq.gz
#> AD26-005-M_S11_MERGED.fastq.gz       AD26-005-M_S11_MERGED.fastq.gz
#> AD30-ABMX-M_S12_MERGED.fastq.gz     AD30-ABMX-M_S12_MERGED.fastq.gz
#> AD32-007-M_S13_MERGED.fastq.gz       AD32-007-M_S13_MERGED.fastq.gz
#> ADABM30X-B_S14_MERGED.fastq.gz       ADABM30X-B_S14_MERGED.fastq.gz
#> ADABM30X-H_S15_MERGED.fastq.gz       ADABM30X-H_S15_MERGED.fastq.gz
#> ADABM30X-M_S16_MERGED.fastq.gz       ADABM30X-M_S16_MERGED.fastq.gz
#> B18-006-B_S19_MERGED.fastq.gz         B18-006-B_S19_MERGED.fastq.gz
#> BA17-050-B_S21_MERGED.fastq.gz       BA17-050-B_S21_MERGED.fastq.gz
#> BB19-006-H_S22_MERGED.fastq.gz       BB19-006-H_S22_MERGED.fastq.gz
#> BB6-019-B_S23_MERGED.fastq.gz         BB6-019-B_S23_MERGED.fastq.gz
#> BB6-019-H_S24_MERGED.fastq.gz         BB6-019-H_S24_MERGED.fastq.gz
#> BB6-019-M_S25_MERGED.fastq.gz         BB6-019-M_S25_MERGED.fastq.gz
#> BE9-006-B_S27_MERGED.fastq.gz         BE9-006-B_S27_MERGED.fastq.gz
#> BE9-006-H_S28_MERGED.fastq.gz         BE9-006-H_S28_MERGED.fastq.gz
#> BE9-006-M_S29_MERGED.fastq.gz         BE9-006-M_S29_MERGED.fastq.gz
#> BG7-010-B_S30_MERGED.fastq.gz         BG7-010-B_S30_MERGED.fastq.gz
#> BG7-010-H_S31_MERGED.fastq.gz         BG7-010-H_S31_MERGED.fastq.gz
#> BG7-010-M_S32_MERGED.fastq.gz         BG7-010-M_S32_MERGED.fastq.gz
#> BJ17-007-M_S34_MERGED.fastq.gz       BJ17-007-M_S34_MERGED.fastq.gz
#> BL7-006-B_S36_MERGED.fastq.gz         BL7-006-B_S36_MERGED.fastq.gz
#> BL7-006-H_S37_MERGED.fastq.gz         BL7-006-H_S37_MERGED.fastq.gz
#> BL7-006-M_S38_MERGED.fastq.gz         BL7-006-M_S38_MERGED.fastq.gz
#> BP11-001-B_S43_MERGED.fastq.gz       BP11-001-B_S43_MERGED.fastq.gz
#> BP11-001-H_S44_MERGED.fastq.gz       BP11-001-H_S44_MERGED.fastq.gz
#> BP11-001-M_S45_MERGED.fastq.gz       BP11-001-M_S45_MERGED.fastq.gz
#> BP12-025-B_S46_MERGED.fastq.gz       BP12-025-B_S46_MERGED.fastq.gz
#> BQ4-018-B_S49_MERGED.fastq.gz         BQ4-018-B_S49_MERGED.fastq.gz
#> BQ4-018-H_S50_MERGED.fastq.gz         BQ4-018-H_S50_MERGED.fastq.gz
#> BQ4-018-M_S51_MERGED.fastq.gz         BQ4-018-M_S51_MERGED.fastq.gz
#> BT-006-M_S55_MERGED.fastq.gz           BT-006-M_S55_MERGED.fastq.gz
#> BV11-002-B_S57_MERGED.fastq.gz       BV11-002-B_S57_MERGED.fastq.gz
#> BV11-002-H_S58_MERGED.fastq.gz       BV11-002-H_S58_MERGED.fastq.gz
#> BV11-002-M_S59_MERGED.fastq.gz       BV11-002-M_S59_MERGED.fastq.gz
#> C21-NV1-B_S62_MERGED.fastq.gz         C21-NV1-B_S62_MERGED.fastq.gz
#> C21-NV1-H_S63_MERGED.fastq.gz         C21-NV1-H_S63_MERGED.fastq.gz
#> C21-NV1-M_S64_MERGED.fastq.gz         C21-NV1-M_S64_MERGED.fastq.gz
#> CB8-019-B_S69_MERGED.fastq.gz         CB8-019-B_S69_MERGED.fastq.gz
#> CB8-019-H_S70_MERGED.fastq.gz         CB8-019-H_S70_MERGED.fastq.gz
#> CB8-019-M_S71_MERGED.fastq.gz         CB8-019-M_S71_MERGED.fastq.gz
#> D18-003-B_S78_MERGED.fastq.gz         D18-003-B_S78_MERGED.fastq.gz
#> D18-003-H_S79_MERGED.fastq.gz         D18-003-H_S79_MERGED.fastq.gz
#> D18-003-M_S80_MERGED.fastq.gz         D18-003-M_S80_MERGED.fastq.gz
#> D61-010-B_S82_MERGED.fastq.gz         D61-010-B_S82_MERGED.fastq.gz
#> D9-027-B_S83_MERGED.fastq.gz           D9-027-B_S83_MERGED.fastq.gz
#> D9-027-H_S84_MERGED.fastq.gz           D9-027-H_S84_MERGED.fastq.gz
#> D9-027-M_S85_MERGED.fastq.gz           D9-027-M_S85_MERGED.fastq.gz
#> DJ2-008-B_S87_MERGED.fastq.gz         DJ2-008-B_S87_MERGED.fastq.gz
#> DJ2-008-H_S88_MERGED.fastq.gz         DJ2-008-H_S88_MERGED.fastq.gz
#> DJ2-008-M_S89_MERGED.fastq.gz         DJ2-008-M_S89_MERGED.fastq.gz
#> DS1-ABM002-B_S91_MERGED.fastq.gz   DS1-ABM002-B_S91_MERGED.fastq.gz
#> DS1-ABM002-H_S92_MERGED.fastq.gz   DS1-ABM002-H_S92_MERGED.fastq.gz
#> DS1-ABM002-M_S93_MERGED.fastq.gz   DS1-ABM002-M_S93_MERGED.fastq.gz
#> DU3-045-B_S94_MERGED.fastq.gz         DU3-045-B_S94_MERGED.fastq.gz
#> DY5-004-B_S96_MERGED.fastq.gz         DY5-004-B_S96_MERGED.fastq.gz
#> DY5-004-H_S97_MERGED.fastq.gz         DY5-004-H_S97_MERGED.fastq.gz
#> DY5-004-M_S98_MERGED.fastq.gz         DY5-004-M_S98_MERGED.fastq.gz
#> E9-009-B_S100_MERGED.fastq.gz         E9-009-B_S100_MERGED.fastq.gz
#> E9-009-H_S101_MERGED.fastq.gz         E9-009-H_S101_MERGED.fastq.gz
#> E9-009-M_S102_MERGED.fastq.gz         E9-009-M_S102_MERGED.fastq.gz
#> EC2-013-B_S104_MERGED.fastq.gz       EC2-013-B_S104_MERGED.fastq.gz
#> F7-015-M_S106_MERGED.fastq.gz         F7-015-M_S106_MERGED.fastq.gz
#> FOMES19-H_S108_MERGED.fastq.gz       FOMES19-H_S108_MERGED.fastq.gz
#> FOMES19-M_S109_MERGED.fastq.gz       FOMES19-M_S109_MERGED.fastq.gz
#> H10-018-M_S110_MERGED.fastq.gz       H10-018-M_S110_MERGED.fastq.gz
#> H24-NVABM1-H_S111_MERGED.fastq.gz H24-NVABM1-H_S111_MERGED.fastq.gz
#> J18-004-B_S114_MERGED.fastq.gz       J18-004-B_S114_MERGED.fastq.gz
#> J18-004-H_S115_MERGED.fastq.gz       J18-004-H_S115_MERGED.fastq.gz
#> J18-004-M_S116_MERGED.fastq.gz       J18-004-M_S116_MERGED.fastq.gz
#> K18-002-H_S117_MERGED.fastq.gz       K18-002-H_S117_MERGED.fastq.gz
#> L19X-B_S119_MERGED.fastq.gz             L19X-B_S119_MERGED.fastq.gz
#> L19X-H_S120_MERGED.fastq.gz             L19X-H_S120_MERGED.fastq.gz
#> L19X-M_S121_MERGED.fastq.gz             L19X-M_S121_MERGED.fastq.gz
#> L23-002-B_S122_MERGED.fastq.gz       L23-002-B_S122_MERGED.fastq.gz
#> L23-002-H_S123_MERGED.fastq.gz       L23-002-H_S123_MERGED.fastq.gz
#> L23-002-M_S124_MERGED.fastq.gz       L23-002-M_S124_MERGED.fastq.gz
#> N19X-B_S126_MERGED.fastq.gz             N19X-B_S126_MERGED.fastq.gz
#> N19X-H_S127_MERGED.fastq.gz             N19X-H_S127_MERGED.fastq.gz
#> N19X-M_S128_MERGED.fastq.gz             N19X-M_S128_MERGED.fastq.gz
#> N22-001-B_S129_MERGED.fastq.gz       N22-001-B_S129_MERGED.fastq.gz
#> N23-002-B_S130_MERGED.fastq.gz       N23-002-B_S130_MERGED.fastq.gz
#> N23-002-H_S131_MERGED.fastq.gz       N23-002-H_S131_MERGED.fastq.gz
#> N23-002-M_S132_MERGED.fastq.gz       N23-002-M_S132_MERGED.fastq.gz
#> NVABM-0163-H_S135_MERGED.fastq.gz NVABM-0163-H_S135_MERGED.fastq.gz
#> NVABM0244-M_S137_MERGED.fastq.gz   NVABM0244-M_S137_MERGED.fastq.gz
#> O20-X-B_S139_MERGED.fastq.gz           O20-X-B_S139_MERGED.fastq.gz
#> O20-X-H_S140_MERGED.fastq.gz           O20-X-H_S140_MERGED.fastq.gz
#> O20-X-M_S141_MERGED.fastq.gz           O20-X-M_S141_MERGED.fastq.gz
#> O21-007-B_S142_MERGED.fastq.gz       O21-007-B_S142_MERGED.fastq.gz
#> O21-007-H_S143_MERGED.fastq.gz       O21-007-H_S143_MERGED.fastq.gz
#> O21-007-M_S144_MERGED.fastq.gz       O21-007-M_S144_MERGED.fastq.gz
#> O24-003-B_S145_MERGED.fastq.gz       O24-003-B_S145_MERGED.fastq.gz
#> O24-003-H_S146_MERGED.fastq.gz       O24-003-H_S146_MERGED.fastq.gz
#> O24-003-M_S147_MERGED.fastq.gz       O24-003-M_S147_MERGED.fastq.gz
#> O26-004-B_S148_MERGED.fastq.gz       O26-004-B_S148_MERGED.fastq.gz
#> O26-004-H_S149_MERGED.fastq.gz       O26-004-H_S149_MERGED.fastq.gz
#> O26-004-M_S150_MERGED.fastq.gz       O26-004-M_S150_MERGED.fastq.gz
#> O9-005-B_S152_MERGED.fastq.gz         O9-005-B_S152_MERGED.fastq.gz
#> P19-023-M_S153_MERGED.fastq.gz       P19-023-M_S153_MERGED.fastq.gz
#> P27-015-M_S154_MERGED.fastq.gz       P27-015-M_S154_MERGED.fastq.gz
#> Q27-ABM003-B_S156_MERGED.fastq.gz Q27-ABM003-B_S156_MERGED.fastq.gz
#> R28-008-B_S158_MERGED.fastq.gz       R28-008-B_S158_MERGED.fastq.gz
#> R28-008-H_S159_MERGED.fastq.gz       R28-008-H_S159_MERGED.fastq.gz
#> R28-008-M_S160_MERGED.fastq.gz       R28-008-M_S160_MERGED.fastq.gz
#> T28-ABM602-B_S162_MERGED.fastq.gz T28-ABM602-B_S162_MERGED.fastq.gz
#> W26-001-B_S165_MERGED.fastq.gz       W26-001-B_S165_MERGED.fastq.gz
#> W26-001-H_S166_MERGED.fastq.gz       W26-001-H_S166_MERGED.fastq.gz
#> W26-001-M_S167_MERGED.fastq.gz       W26-001-M_S167_MERGED.fastq.gz
#> W9-025-M_S169_MERGED.fastq.gz         W9-025-M_S169_MERGED.fastq.gz
#> X24-009-B_S170_MERGED.fastq.gz       X24-009-B_S170_MERGED.fastq.gz
#> X24-009-H_S171_MERGED.fastq.gz       X24-009-H_S171_MERGED.fastq.gz
#> X24-009-M_S172_MERGED.fastq.gz       X24-009-M_S172_MERGED.fastq.gz
#> X29-004-B_S174_MERGED.fastq.gz       X29-004-B_S174_MERGED.fastq.gz
#> X29-004-H_S175_MERGED.fastq.gz       X29-004-H_S175_MERGED.fastq.gz
#> X29-004-M_S176_MERGED.fastq.gz       X29-004-M_S176_MERGED.fastq.gz
#> Y21-ABM484-H_S177_MERGED.fastq.gz Y21-ABM484-H_S177_MERGED.fastq.gz
#> Y28-002-B_S178_MERGED.fastq.gz       Y28-002-B_S178_MERGED.fastq.gz
#> Y28-002-H_S179_MERGED.fastq.gz       Y28-002-H_S179_MERGED.fastq.gz
#> Y28-002-M_S180_MERGED.fastq.gz       Y28-002-M_S180_MERGED.fastq.gz
#> Y29-007-B_S181_MERGED.fastq.gz       Y29-007-B_S181_MERGED.fastq.gz
#> Y29-007-H_S182_MERGED.fastq.gz       Y29-007-H_S182_MERGED.fastq.gz
#> Y29-007-M_S183_MERGED.fastq.gz       Y29-007-M_S183_MERGED.fastq.gz
#> Y31-ABM484-B_S184_MERGED.fastq.gz Y31-ABM484-B_S184_MERGED.fastq.gz
#> Z29-001-H_S185_MERGED.fastq.gz       Z29-001-H_S185_MERGED.fastq.gz
#> Z30-ABM560-M_S187_MERGED.fastq.gz Z30-ABM560-M_S187_MERGED.fastq.gz
#> 
#> $dabest_objects
#> $dabest_objects$Hill_0
#> DABESTR v2025.3.15
#> ==================
#> 
#> Good afternoon!
#> The current time is 17:19 PM on Monday June 15, 2026.
#> 
#> The character(0) Cohen's d between Low and High is -0.058 [95%CI -0.487, 0.371].
#> The p-value of the two-sided permutation t-test is 0.7902, calculated for legacy purposes only.
#> 
#> The character(0) Cohen's d between Middle and High is -0.221 [95%CI -0.622, 0.211].
#> The p-value of the two-sided permutation t-test is 0.3133, calculated for legacy purposes only.
#> 
#> 5000 bootstrap samples were taken; the confidence interval is bias-corrected and accelerated.
#> Any p-value reported is the probability of observing the effect size (or greater),
#> assuming the null hypothesis of zero difference is true.
#> For each p-value, 5000 reshuffles of the control and test labels were performed.
#> 
#> 
#> $dabest_objects$Hill_1
#> DABESTR v2025.3.15
#> ==================
#> 
#> Good afternoon!
#> The current time is 17:19 PM on Monday June 15, 2026.
#> 
#> The character(0) Cohen's d between Low and High is 0.012 [95%CI -0.412, 0.429].
#> The p-value of the two-sided permutation t-test is 0.9566, calculated for legacy purposes only.
#> 
#> The character(0) Cohen's d between Middle and High is 0.105 [95%CI -0.316, 0.529].
#> The p-value of the two-sided permutation t-test is 0.6248, calculated for legacy purposes only.
#> 
#> 5000 bootstrap samples were taken; the confidence interval is bias-corrected and accelerated.
#> Any p-value reported is the probability of observing the effect size (or greater),
#> assuming the null hypothesis of zero difference is true.
#> For each p-value, 5000 reshuffles of the control and test labels were performed.
#> 
#> 
#> $dabest_objects$Hill_2
#> DABESTR v2025.3.15
#> ==================
#> 
#> Good afternoon!
#> The current time is 17:19 PM on Monday June 15, 2026.
#> 
#> The character(0) Cohen's d between Low and High is -0.033 [95%CI -0.467, 0.393].
#> The p-value of the two-sided permutation t-test is 0.8770, calculated for legacy purposes only.
#> 
#> The character(0) Cohen's d between Middle and High is 0.128 [95%CI -0.291, 0.541].
#> The p-value of the two-sided permutation t-test is 0.5518, calculated for legacy purposes only.
#> 
#> 5000 bootstrap samples were taken; the confidence interval is bias-corrected and accelerated.
#> Any p-value reported is the probability of observing the effect size (or greater),
#> assuming the null hypothesis of zero difference is true.
#> For each p-value, 5000 reshuffles of the control and test labels were performed.
#> 
#> 
#> 
#> $plots
#> $plots$Hill_0

#> 
#> $plots$Hill_1

#> 
#> $plots$Hill_2

#> 
#> 
#> $summary
#> # A tibble: 6 × 8
#>   metric comparison   effect_size ci_lower ci_upper pvalue_permtest pvalue_welch
#>   <chr>  <chr>              <dbl>    <dbl>    <dbl>           <dbl>        <dbl>
#> 1 Hill_0 High vs Low      -0.0581   -0.487    0.371           0.780        0.790
#> 2 Hill_0 High vs Mid…     -0.221    -0.622    0.211           0.310        0.313
#> 3 Hill_1 High vs Low       0.0118   -0.412    0.429           0.959        0.957
#> 4 Hill_1 High vs Mid…      0.105    -0.316    0.529           0.624        0.625
#> 5 Hill_2 High vs Low      -0.0334   -0.467    0.393           0.88         0.877
#> 6 Hill_2 High vs Mid…      0.128    -0.291    0.541           0.557        0.552
#> # ℹ 1 more variable: pvalue_mann_whitney <dbl>
#> 
#> $effect_type
#> [1] "cohens_d"
#> 
#> attr(,"class")
#> [1] "estim_diff_pq_result"
res$plots$Hill_0

res$summary
#> # A tibble: 6 × 8
#>   metric comparison   effect_size ci_lower ci_upper pvalue_permtest pvalue_welch
#>   <chr>  <chr>              <dbl>    <dbl>    <dbl>           <dbl>        <dbl>
#> 1 Hill_0 High vs Low      -0.0581   -0.487    0.371           0.780        0.790
#> 2 Hill_0 High vs Mid…     -0.221    -0.622    0.211           0.310        0.313
#> 3 Hill_1 High vs Low       0.0118   -0.412    0.429           0.959        0.957
#> 4 Hill_1 High vs Mid…      0.105    -0.316    0.529           0.624        0.625
#> 5 Hill_2 High vs Low      -0.0334   -0.467    0.393           0.88         0.877
#> 6 Hill_2 High vs Mid…      0.128    -0.291    0.541           0.557        0.552
#> # ℹ 1 more variable: pvalue_mann_whitney <dbl>
```

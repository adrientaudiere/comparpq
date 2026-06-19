# ANCOM-BC analysis on each phyloseq object in a list_phyloseq

Performs ANCOM-BC differential abundance analysis using
[`MiscMetabar::ancombc_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ancombc_pq.html)
on each phyloseq object in a list_phyloseq and returns a combined result
table.

## Usage

``` r
ancombc_lpq(x, fact, levels_fact = NULL, tax_level = NULL, verbose = TRUE, ...)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- fact:

  (character, required) The name of a column in `sample_data` to use as
  the grouping factor. Must be present in **all** phyloseq objects.

- levels_fact:

  (character vector, default NULL) Levels of the factor to include. If
  NULL, all levels are used.

- tax_level:

  (character, default "Class") Taxonomic level for agglomeration before
  analysis.

- verbose:

  (logical, default TRUE) If TRUE, print progress messages.

- ...:

  Additional arguments passed to
  [`MiscMetabar::ancombc_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ancombc_pq.html).

## Value

A tibble with the combined ANCOM-BC results from all phyloseq objects.
Contains the `$res` data frame from each ANCOM-BC run, with an
additional `name` column identifying the source phyloseq object.

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This function requires that the list_phyloseq type is NOT
`SEPARATE_ANALYSIS`, as the factor must be common across all phyloseq
objects.

When no common taxa exist across the phyloseq objects, taxa names are
suffixed with the phyloseq object name to make them distinguishable.

## See also

[`MiscMetabar::ancombc_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ancombc_pq.html),
[`aldex_lpq()`](https://adrientaudiere.github.io/comparpq/reference/aldex_lpq.md),
[`multipatt_lpq()`](https://adrientaudiere.github.io/comparpq/reference/multipatt_lpq.md)

## Examples

``` r
# Subset to the 80 most abundant taxa to keep the example fast
# (the full data_fungi has 1420 taxa, which is slow for ANCOMBC).
data_fungi_small <- prune_taxa(
  names(sort(taxa_sums(data_fungi), decreasing = TRUE))[1:80],
  data_fungi
)
data_fungi_small <- clean_pq(prune_samples(
  sample_sums(data_fungi_small) >= 500, data_fungi_small
))
data_fungi_high <- multiply_counts_pq(data_fungi_small, "Height", "High", 2)
#> Modified 40 taxa in 28 matched samples

lpq <- list_phyloseq(
  list(
    fungi = data_fungi_small,
    fungi_height = data_fungi_high
  )
)
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: REPRODUCIBILITY
#> ℹ 127 common samples, 80 common taxa
#> ✔ list_phyloseq created (REPRODUCIBILITY)

results <- ancombc_lpq(lpq, fact = "Height")
#> Running ANCOM-BC on 2 phyloseq objects
#> Factor: Height | Tax level: 
#> → Processing: fungi
#> Checking the input data type ...
#> The input data is of type: TreeSummarizedExperiment
#> PASS
#> Checking the sample metadata ...
#> The specified variables in the formula: Height
#> The available variables in the sample metadata: X, Sample_names, Tree_name, Sample_id, Height, Diameter, Time
#> PASS
#> Checking other arguments ...
#> The number of groups of interest is: 3
#> The sample size per group is: High = 28, Low = 26, Middle = 27
#> PASS
#> Warning: The number of taxa used for estimating sample-specific biases is: 33
#> A large number of taxa (>50) is required for the consistent estimation of biases
#> Obtaining initial estimates ...
#> Estimating sample-specific biases ...
#> Loading required package: foreach
#> Loading required package: rngtools
#> ANCOM-BC2 primary results ...
#> Conducting sensitivity analysis for pseudo-count addition to 0s ...
#> For taxa that are significant but do not pass the sensitivity analysis,
#> they are marked in the 'passed_ss' column and will be treated as non-significant in the 'diff_robust' column.
#> For detailed instructions on performing sensitivity analysis, please refer to the package vignette.
#> → Processing: fungi_height
#> Checking the input data type ...
#> The input data is of type: TreeSummarizedExperiment
#> PASS
#> Checking the sample metadata ...
#> The specified variables in the formula: Height
#> The available variables in the sample metadata: X, Sample_names, Tree_name, Sample_id, Height, Diameter, Time
#> PASS
#> Checking other arguments ...
#> The number of groups of interest is: 3
#> The sample size per group is: High = 28, Low = 26, Middle = 27
#> PASS
#> Warning: The number of taxa used for estimating sample-specific biases is: 33
#> A large number of taxa (>50) is required for the consistent estimation of biases
#> Obtaining initial estimates ...
#> Estimating sample-specific biases ...
#> ANCOM-BC2 primary results ...
#> Conducting sensitivity analysis for pseudo-count addition to 0s ...
#> For taxa that are significant but do not pass the sensitivity analysis,
#> they are marked in the 'passed_ss' column and will be treated as non-significant in the 'diff_robust' column.
#> For detailed instructions on performing sensitivity analysis, please refer to the package vignette.
results
#> # A tibble: 66 × 26
#>    name  taxon `lfc_(Intercept)` lfc_HeightLow lfc_HeightMiddle `se_(Intercept)`
#>    <chr> <chr>             <dbl>         <dbl>            <dbl>            <dbl>
#>  1 fungi ASV2             0.596        -0.105            -0.185            0.686
#>  2 fungi ASV6             1.56         -1.69             -1.12             0.518
#>  3 fungi ASV7             0.633        -0.709            -0.383            0.529
#>  4 fungi ASV8             0.0438       -0.148             1.87             0.728
#>  5 fungi ASV10            0.462        -0.936            -0.164            0.539
#>  6 fungi ASV12            0.0906       -0.0560           -0.240            0.647
#>  7 fungi ASV13           -0.544         0.497             0.708            0.601
#>  8 fungi ASV18           -0.108        -0.972            -0.200            0.524
#>  9 fungi ASV19            0.0388        1.46              1.28             0.470
#> 10 fungi ASV24           -2.89          3.76              0.458            0.516
#> # ℹ 56 more rows
#> # ℹ 20 more variables: se_HeightLow <dbl>, se_HeightMiddle <dbl>,
#> #   `W_(Intercept)` <dbl>, W_HeightLow <dbl>, W_HeightMiddle <dbl>,
#> #   `p_(Intercept)` <dbl>, p_HeightLow <dbl>, p_HeightMiddle <dbl>,
#> #   `q_(Intercept)` <dbl>, q_HeightLow <dbl>, q_HeightMiddle <dbl>,
#> #   `diff_(Intercept)` <lgl>, diff_HeightLow <lgl>, diff_HeightMiddle <lgl>,
#> #   `passed_ss_(Intercept)` <lgl>, passed_ss_HeightLow <lgl>, …

results |>
  filter(diff_HeightLow) |>
  ggplot(aes(y = taxon, x = lfc_HeightLow, color = name, shape = q_HeightLow < 0.1)) +
  geom_point(size = 4) +
  geom_vline(xintercept = 0) +
  theme_minimal()


results_Genus <- ancombc_lpq(lpq, fact = "Height", tax_level = "Genus")
#> Running ANCOM-BC on 2 phyloseq objects
#> Factor: Height | Tax level: Genus
#> → Processing: fungi
#> Checking the input data type ...
#> The input data is of type: TreeSummarizedExperiment
#> PASS
#> Checking the sample metadata ...
#> The specified variables in the formula: Height
#> The available variables in the sample metadata: X, Sample_names, Tree_name, Sample_id, Height, Diameter, Time
#> PASS
#> Checking other arguments ...
#> The number of groups of interest is: 3
#> The sample size per group is: High = 28, Low = 26, Middle = 27
#> PASS
#> Warning: The number of taxa used for estimating sample-specific biases is: 33
#> A large number of taxa (>50) is required for the consistent estimation of biases
#> Obtaining initial estimates ...
#> Estimating sample-specific biases ...
#> ANCOM-BC2 primary results ...
#> Conducting sensitivity analysis for pseudo-count addition to 0s ...
#> For taxa that are significant but do not pass the sensitivity analysis,
#> they are marked in the 'passed_ss' column and will be treated as non-significant in the 'diff_robust' column.
#> For detailed instructions on performing sensitivity analysis, please refer to the package vignette.
#> → Processing: fungi_height
#> Checking the input data type ...
#> The input data is of type: TreeSummarizedExperiment
#> PASS
#> Checking the sample metadata ...
#> The specified variables in the formula: Height
#> The available variables in the sample metadata: X, Sample_names, Tree_name, Sample_id, Height, Diameter, Time
#> PASS
#> Checking other arguments ...
#> The number of groups of interest is: 3
#> The sample size per group is: High = 28, Low = 26, Middle = 27
#> PASS
#> Warning: The number of taxa used for estimating sample-specific biases is: 33
#> A large number of taxa (>50) is required for the consistent estimation of biases
#> Obtaining initial estimates ...
#> Estimating sample-specific biases ...
#> ANCOM-BC2 primary results ...
#> Conducting sensitivity analysis for pseudo-count addition to 0s ...
#> For taxa that are significant but do not pass the sensitivity analysis,
#> they are marked in the 'passed_ss' column and will be treated as non-significant in the 'diff_robust' column.
#> For detailed instructions on performing sensitivity analysis, please refer to the package vignette.

results_Genus |>
  filter(diff_HeightLow) |>
  ggplot(aes(y = taxon, x = lfc_HeightLow, color = name, shape = q_HeightLow < 0.1)) +
  geom_point(size = 4) +
  geom_vline(xintercept = 0) +
  theme_minimal()

```

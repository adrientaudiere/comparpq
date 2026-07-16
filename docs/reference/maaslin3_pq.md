# Run MaAsLin3 differential abundance analysis on a phyloseq object

A wrapper around
[`maaslin3::maaslin3()`](https://rdrr.io/pkg/maaslin3/man/maaslin3.html)
for phyloseq objects. Optionally includes the number of reads (library
size) as a covariate in the model to account for differences in
sequencing depth.

## Usage

``` r
maaslin3_pq(
  physeq,
  formula,
  reference = NULL,
  correction_for_sample_size = TRUE,
  output = "res_maaslin3",
  ...
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- formula:

  (character, required) A formula string for the model (e.g.,
  `"~ Height"` or `"~ Height + Time"`).

- reference:

  (named list, default NULL) Reference levels for categorical variables
  with more than 2 levels. Format: `list(varname = "ref_level")`. For
  example: `list(Height = "Low")` sets "Low" as the reference for
  Height. Variables not in this list will use their first level as
  reference.

- correction_for_sample_size:

  (logical, default TRUE) If TRUE, adds `nb_seq` (library size) to the
  formula to control for sequencing depth.

- output:

  (character, default "res_maaslin3") Output directory for MaAsLin3
  results.

- ...:

  Additional arguments passed to
  [`maaslin3::maaslin3()`](https://rdrr.io/pkg/maaslin3/man/maaslin3.html).

## Value

The result object from
[`maaslin3::maaslin3()`](https://rdrr.io/pkg/maaslin3/man/maaslin3.html),
containing:

- `results`: Data frame with differential abundance results

- `results_ordered`: Results ordered by significance

- Other MaAsLin3 output components

## Details

\#' \#TODO VERY experimental
[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

MaAsLin3 requires categorical variables with more than 2 levels to have
a defined reference level. This function handles this by:

1.  Converting character variables to factors

2.  Setting reference levels via the `reference` parameter

3.  Using [`relevel()`](https://rdrr.io/r/stats/relevel.html) to set the
    specified reference as the first level

The `correction_for_sample_size` option adds library size as a
covariate, which can help account for compositional effects and
sequencing depth differences between samples.

## References

Nickols WA, et al. (2024). MaAsLin 3: Refining and extending generalized
multivariable linear models for meta-omic association discovery.
bioRxiv.

## See also

[`MiscMetabar::aldex_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/aldex_pq.html),
[`ancombc_lpq()`](https://adrientaudiere.github.io/comparpq/reference/ancombc_lpq.md),
[`gg_maaslin3_plot()`](https://adrientaudiere.github.io/comparpq/reference/gg_maaslin3_plot.md)

## Author

Adrien Taudière

## Examples

``` r
# Basic usage with a binary factor (45 taxa, fast)
res <- maaslin3_pq(data_fungi_mini, formula = "~ Height")
#> Taxa are now in columns.
#> 2026-07-16 12:41:00.98 INFO::Writing function arguments to log file
#> 2026-07-16 12:41:01.04 INFO::Verifying options selected are valid
#> 2026-07-16 12:41:01.04 INFO::Determining format of input files
#> 2026-07-16 12:41:01.04 INFO::Input format is data samples as rows and metadata samples as rows
#> 2026-07-16 12:41:01.04 INFO::Running selected normalization method: TSS
#> 2026-07-16 12:41:01.04 INFO::Writing normalized data to file res_maaslin3/features/data_norm.tsv
#> 2026-07-16 12:41:01.05 INFO::Filter data based on min abundance, min prevalence, and max prevalence
#> 2026-07-16 12:41:01.05 INFO::Total samples in data: 137
#> 2026-07-16 12:41:01.05 INFO::Min samples required with min abundance for a feature not to be filtered: 0.000000
#> 2026-07-16 12:41:01.05 INFO::Max samples allowed with min abundance for a feature not to be filtered: 138.370000
#> 2026-07-16 12:41:01.06 INFO::Total filtered features: 0
#> 2026-07-16 12:41:01.06 INFO::Filtered feature names from abundance, min prevalence, and max prevalence filtering:
#> 2026-07-16 12:41:01.06 INFO::Total features filtered by non-zero variance filtering: 2
#> 2026-07-16 12:41:01.06 INFO::Filtered feature names from variance filtering: ASV54, ASV108
#> 2026-07-16 12:41:01.06 INFO::Writing filtered data to file res_maaslin3/features/filtered_data.tsv
#> 2026-07-16 12:41:01.06 INFO::Running selected transform method: LOG
#> 2026-07-16 12:41:01.07 INFO::Writing normalized, filtered, transformed data to file res_maaslin3/features/data_transformed.tsv
#> 2026-07-16 12:41:01.07 INFO::Factor detected for categorial metadata 'Height'. Using as-is.
#> 2026-07-16 12:41:01.07 INFO::Applying z-score to standardize continuous metadata
#> 2026-07-16 12:41:01.10 INFO::Running the linear model component
#> 2026-07-16 12:41:01.12 INFO::Fitting model to feature number 1, ASV7
#> 2026-07-16 12:41:01.12 INFO::Fitting model to feature number 2, ASV8
#> 2026-07-16 12:41:01.13 INFO::Fitting model to feature number 3, ASV12
#> 2026-07-16 12:41:01.13 INFO::Fitting model to feature number 4, ASV18
#> 2026-07-16 12:41:01.14 INFO::Fitting model to feature number 5, ASV25
#> 2026-07-16 12:41:01.14 INFO::Fitting model to feature number 6, ASV26
#> 2026-07-16 12:41:01.15 INFO::Fitting model to feature number 7, ASV27
#> 2026-07-16 12:41:01.15 INFO::Fitting model to feature number 8, ASV29
#> 2026-07-16 12:41:01.16 INFO::Fitting model to feature number 9, ASV32
#> 2026-07-16 12:41:01.16 INFO::Fitting model to feature number 10, ASV34
#> 2026-07-16 12:41:01.16 INFO::Fitting model to feature number 11, ASV35
#> 2026-07-16 12:41:01.17 INFO::Fitting model to feature number 12, ASV41
#> 2026-07-16 12:41:01.17 INFO::Fitting model to feature number 13, ASV42
#> 2026-07-16 12:41:01.17 INFO::Fitting model to feature number 14, ASV46
#> 2026-07-16 12:41:01.18 INFO::Fitting model to feature number 15, ASV47
#> 2026-07-16 12:41:01.18 INFO::Fitting model to feature number 16, ASV48
#> 2026-07-16 12:41:01.18 WARNING::Fitting problem for feature 16 returning NA
#> 2026-07-16 12:41:01.18 INFO::Fitting model to feature number 17, ASV49
#> 2026-07-16 12:41:01.19 INFO::Fitting model to feature number 18, ASV50
#> 2026-07-16 12:41:01.19 WARNING::Fitting problem for feature 18 returning NA
#> 2026-07-16 12:41:01.19 INFO::Fitting model to feature number 19, ASV53
#> 2026-07-16 12:41:01.20 INFO::Fitting model to feature number 20, ASV58
#> 2026-07-16 12:41:01.20 INFO::Fitting model to feature number 21, ASV59
#> 2026-07-16 12:41:01.20 INFO::Fitting model to feature number 22, ASV61
#> 2026-07-16 12:41:01.21 INFO::Fitting model to feature number 23, ASV62
#> 2026-07-16 12:41:01.21 INFO::Fitting model to feature number 24, ASV63
#> 2026-07-16 12:41:01.22 INFO::Fitting model to feature number 25, ASV64
#> 2026-07-16 12:41:01.22 INFO::Fitting model to feature number 26, ASV67
#> 2026-07-16 12:41:01.22 INFO::Fitting model to feature number 27, ASV68
#> 2026-07-16 12:41:01.22 INFO::Fitting model to feature number 28, ASV71
#> 2026-07-16 12:41:01.22 INFO::Fitting model to feature number 29, ASV72
#> 2026-07-16 12:41:01.23 INFO::Fitting model to feature number 30, ASV75
#> 2026-07-16 12:41:01.23 INFO::Fitting model to feature number 31, ASV77
#> 2026-07-16 12:41:01.24 WARNING::Fitting problem for feature 31 returning NA
#> 2026-07-16 12:41:01.24 INFO::Fitting model to feature number 32, ASV82
#> 2026-07-16 12:41:01.24 INFO::Fitting model to feature number 33, ASV83
#> 2026-07-16 12:41:01.24 INFO::Fitting model to feature number 34, ASV85
#> 2026-07-16 12:41:01.25 INFO::Fitting model to feature number 35, ASV91
#> 2026-07-16 12:41:01.25 INFO::Fitting model to feature number 36, ASV93
#> 2026-07-16 12:41:01.26 WARNING::Fitting problem for feature 36 returning NA
#> 2026-07-16 12:41:01.26 INFO::Fitting model to feature number 37, ASV94
#> 2026-07-16 12:41:01.26 INFO::Fitting model to feature number 38, ASV99
#> 2026-07-16 12:41:01.27 INFO::Fitting model to feature number 39, ASV100
#> 2026-07-16 12:41:01.27 INFO::Fitting model to feature number 40, ASV101
#> 2026-07-16 12:41:01.27 INFO::Fitting model to feature number 41, ASV104
#> 2026-07-16 12:41:01.28 INFO::Fitting model to feature number 42, ASV105
#> 2026-07-16 12:41:01.28 INFO::Fitting model to feature number 43, ASV107
#> 2026-07-16 12:41:01.29 INFO::Performing tests against medians
#> 2026-07-16 12:41:02.21 INFO::Counting total values for each feature
#> 2026-07-16 12:41:02.22 INFO::Running the logistic model component
#> 2026-07-16 12:41:02.24 INFO::Fitting model to feature number 1, ASV7
#> 2026-07-16 12:41:02.25 INFO::Fitting model to feature number 2, ASV8
#> 2026-07-16 12:41:02.25 INFO::Fitting model to feature number 3, ASV12
#> 2026-07-16 12:41:02.26 INFO::Fitting model to feature number 4, ASV18
#> 2026-07-16 12:41:02.27 INFO::Fitting model to feature number 5, ASV25
#> 2026-07-16 12:41:02.27 INFO::Fitting model to feature number 6, ASV26
#> 2026-07-16 12:41:02.28 INFO::Fitting model to feature number 7, ASV27
#> 2026-07-16 12:41:02.29 INFO::Fitting model to feature number 8, ASV29
#> 2026-07-16 12:41:02.30 INFO::Fitting model to feature number 9, ASV32
#> 2026-07-16 12:41:02.30 INFO::Fitting model to feature number 10, ASV34
#> 2026-07-16 12:41:02.31 INFO::Fitting model to feature number 11, ASV35
#> 2026-07-16 12:41:02.31 INFO::Fitting model to feature number 12, ASV41
#> 2026-07-16 12:41:02.32 INFO::Fitting model to feature number 13, ASV42
#> 2026-07-16 12:41:02.33 INFO::Fitting model to feature number 14, ASV46
#> 2026-07-16 12:41:02.33 INFO::Fitting model to feature number 15, ASV47
#> 2026-07-16 12:41:02.34 INFO::Fitting model to feature number 16, ASV48
#> 2026-07-16 12:41:02.35 INFO::Fitting model to feature number 17, ASV49
#> 2026-07-16 12:41:02.35 INFO::Fitting model to feature number 18, ASV50
#> 2026-07-16 12:41:02.36 INFO::Fitting model to feature number 19, ASV53
#> 2026-07-16 12:41:02.36 INFO::Fitting model to feature number 20, ASV58
#> 2026-07-16 12:41:02.37 INFO::Fitting model to feature number 21, ASV59
#> 2026-07-16 12:41:02.37 INFO::Fitting model to feature number 22, ASV61
#> 2026-07-16 12:41:02.38 INFO::Fitting model to feature number 23, ASV62
#> 2026-07-16 12:41:02.39 INFO::Fitting model to feature number 24, ASV63
#> 2026-07-16 12:41:02.40 INFO::Fitting model to feature number 25, ASV64
#> 2026-07-16 12:41:02.40 INFO::Fitting model to feature number 26, ASV67
#> 2026-07-16 12:41:02.41 INFO::Fitting model to feature number 27, ASV68
#> 2026-07-16 12:41:02.41 INFO::Fitting model to feature number 28, ASV71
#> 2026-07-16 12:41:02.42 INFO::Fitting model to feature number 29, ASV72
#> 2026-07-16 12:41:02.43 INFO::Fitting model to feature number 30, ASV75
#> 2026-07-16 12:41:02.43 INFO::Fitting model to feature number 31, ASV77
#> 2026-07-16 12:41:02.44 INFO::Fitting model to feature number 32, ASV82
#> 2026-07-16 12:41:02.45 INFO::Fitting model to feature number 33, ASV83
#> 2026-07-16 12:41:02.45 INFO::Fitting model to feature number 34, ASV85
#> 2026-07-16 12:41:02.46 INFO::Fitting model to feature number 35, ASV91
#> 2026-07-16 12:41:02.46 INFO::Fitting model to feature number 36, ASV93
#> 2026-07-16 12:41:02.47 INFO::Fitting model to feature number 37, ASV94
#> 2026-07-16 12:41:02.48 INFO::Fitting model to feature number 38, ASV99
#> 2026-07-16 12:41:02.48 INFO::Fitting model to feature number 39, ASV100
#> 2026-07-16 12:41:02.49 INFO::Fitting model to feature number 40, ASV101
#> 2026-07-16 12:41:02.49 INFO::Fitting model to feature number 41, ASV104
#> 2026-07-16 12:41:02.50 INFO::Fitting model to feature number 42, ASV105
#> 2026-07-16 12:41:02.51 INFO::Fitting model to feature number 43, ASV107
#> 2026-07-16 12:41:02.52 INFO::Counting total values for each feature
#> 2026-07-16 12:41:02.52 INFO::Re-running abundances for warn_prevalence
#> 2026-07-16 12:41:02.52 INFO::Running selected normalization method: TSS
#> 2026-07-16 12:41:02.52 INFO::Running selected transform method: LOG
#> 2026-07-16 12:41:02.54 INFO::Fitting model to feature number 1, ASV7
#> 2026-07-16 12:41:02.55 INFO::Fitting model to feature number 2, ASV8
#> 2026-07-16 12:41:02.55 INFO::Fitting model to feature number 3, ASV12
#> 2026-07-16 12:41:02.56 INFO::Fitting model to feature number 4, ASV18
#> 2026-07-16 12:41:02.56 INFO::Fitting model to feature number 5, ASV25
#> 2026-07-16 12:41:02.57 INFO::Fitting model to feature number 6, ASV26
#> 2026-07-16 12:41:02.57 INFO::Fitting model to feature number 7, ASV27
#> 2026-07-16 12:41:02.57 INFO::Fitting model to feature number 8, ASV29
#> 2026-07-16 12:41:02.58 INFO::Fitting model to feature number 9, ASV32
#> 2026-07-16 12:41:02.58 INFO::Fitting model to feature number 10, ASV34
#> 2026-07-16 12:41:02.59 INFO::Fitting model to feature number 11, ASV35
#> 2026-07-16 12:41:02.59 INFO::Fitting model to feature number 12, ASV41
#> 2026-07-16 12:41:02.60 INFO::Fitting model to feature number 13, ASV42
#> 2026-07-16 12:41:02.60 INFO::Fitting model to feature number 14, ASV46
#> 2026-07-16 12:41:02.60 INFO::Fitting model to feature number 15, ASV47
#> 2026-07-16 12:41:02.61 INFO::Fitting model to feature number 16, ASV48
#> 2026-07-16 12:41:02.61 WARNING::Fitting problem for feature 16 returning NA
#> 2026-07-16 12:41:02.61 INFO::Fitting model to feature number 17, ASV49
#> 2026-07-16 12:41:02.61 INFO::Fitting model to feature number 18, ASV50
#> 2026-07-16 12:41:02.62 WARNING::Fitting problem for feature 18 returning NA
#> 2026-07-16 12:41:02.62 INFO::Fitting model to feature number 19, ASV53
#> 2026-07-16 12:41:02.62 INFO::Fitting model to feature number 20, ASV58
#> 2026-07-16 12:41:02.63 INFO::Fitting model to feature number 21, ASV59
#> 2026-07-16 12:41:02.63 INFO::Fitting model to feature number 22, ASV61
#> 2026-07-16 12:41:02.64 INFO::Fitting model to feature number 23, ASV62
#> 2026-07-16 12:41:02.64 INFO::Fitting model to feature number 24, ASV63
#> 2026-07-16 12:41:02.65 INFO::Fitting model to feature number 25, ASV64
#> 2026-07-16 12:41:02.65 INFO::Fitting model to feature number 26, ASV67
#> 2026-07-16 12:41:02.65 INFO::Fitting model to feature number 27, ASV68
#> 2026-07-16 12:41:02.65 INFO::Fitting model to feature number 28, ASV71
#> 2026-07-16 12:41:02.66 INFO::Fitting model to feature number 29, ASV72
#> 2026-07-16 12:41:02.66 INFO::Fitting model to feature number 30, ASV75
#> 2026-07-16 12:41:02.67 INFO::Fitting model to feature number 31, ASV77
#> 2026-07-16 12:41:02.67 WARNING::Fitting problem for feature 31 returning NA
#> 2026-07-16 12:41:02.67 INFO::Fitting model to feature number 32, ASV82
#> 2026-07-16 12:41:02.67 INFO::Fitting model to feature number 33, ASV83
#> 2026-07-16 12:41:02.68 INFO::Fitting model to feature number 34, ASV85
#> 2026-07-16 12:41:02.68 INFO::Fitting model to feature number 35, ASV91
#> 2026-07-16 12:41:02.69 INFO::Fitting model to feature number 36, ASV93
#> 2026-07-16 12:41:02.69 WARNING::Fitting problem for feature 36 returning NA
#> 2026-07-16 12:41:02.69 INFO::Fitting model to feature number 37, ASV94
#> 2026-07-16 12:41:02.70 INFO::Fitting model to feature number 38, ASV99
#> 2026-07-16 12:41:02.73 INFO::Fitting model to feature number 39, ASV100
#> 2026-07-16 12:41:02.74 INFO::Fitting model to feature number 40, ASV101
#> 2026-07-16 12:41:02.74 INFO::Fitting model to feature number 41, ASV104
#> 2026-07-16 12:41:02.75 INFO::Fitting model to feature number 42, ASV105
#> 2026-07-16 12:41:02.75 INFO::Fitting model to feature number 43, ASV107
#> 2026-07-16 12:41:02.78 WARNING::Deleting existing residuals file: res_maaslin3/fits/residuals_linear.rds
#> 2026-07-16 12:41:02.78 INFO::Writing residuals to file res_maaslin3/fits/residuals_linear.rds
#> 2026-07-16 12:41:02.78 WARNING::Deleting existing fitted file: res_maaslin3/fits/fitted_linear.rds
#> 2026-07-16 12:41:02.78 INFO::Writing fitted values to file res_maaslin3/fits/fitted_linear.rds
#> 2026-07-16 12:41:02.78 WARNING::Deleting existing residuals file: res_maaslin3/fits/residuals_logistic.rds
#> 2026-07-16 12:41:02.78 INFO::Writing residuals to file res_maaslin3/fits/residuals_logistic.rds
#> 2026-07-16 12:41:02.79 WARNING::Deleting existing fitted file: res_maaslin3/fits/fitted_logistic.rds
#> 2026-07-16 12:41:02.79 INFO::Writing fitted values to file res_maaslin3/fits/fitted_logistic.rds
#> 2026-07-16 12:41:02.79 INFO::Writing all the results to file (ordered 
#>             by increasing individual q-values): res_maaslin3/all_results.tsv
#> 2026-07-16 12:41:02.80 INFO::Writing the significant results without errors (those which have joint q-values less than or equal to the threshold of 0.100000 ) to file (ordered by increasing individual q-values): res_maaslin3/significant_results.tsv
#> 2026-07-16 12:41:02.80 INFO::Writing summary plot of significant
#>                         results to file: res_maaslin3/figures/summary_plot.pdf
#> 2026-07-16 12:41:04.90 INFO::Writing association plots (one for each significant association) to output folder: res_maaslin3/figures
#> 2026-07-16 12:41:04.91 INFO::All associations had errors 
#>                                 or were insignificant.
res$results
#> NULL
if (FALSE) { # \dontrun{
# Specify reference level for multi-level factor
res <- maaslin3_pq(
  data_fungi,
  formula = "~ Height",
  reference = list(Height = "Low")
)

# Multiple variables with references
res <- maaslin3_pq(
  data_fungi,
  formula = "~ Height + Site",
  reference = list(Height = "Low", Site = "A")
)

# Without library size correction
res <- maaslin3_pq(
  data_fungi,
  formula = "~ Height",
  correction_for_sample_size = FALSE
)

# Plot results
gg_maaslin3_plot(res, type = "volcano")

# Full example with GlobalPatterns dataset
data("GlobalPatterns")

# Subset to two very different environments: Feces vs Soil
gp_subset <- subset_samples(GlobalPatterns, SampleType %in% c("Feces", "Soil"))

# Agglomerate at Phylum level for clearer results
gp_phylum <- tax_glom(gp_subset, taxrank = "Phylum", NArm = FALSE)

# Run MaAsLin3 with Soil as reference
res <- maaslin3_pq(
  gp_phylum,
  formula = "~ SampleType",
  reference = list(SampleType = "Soil"),
  output = "output/maaslin3_example",
  correction_for_sample_size = FALSE
)

# Plot results
gg_maaslin3_plot(res, type = "volcano", signif_threshold = 0.1)
gg_maaslin3_plot(res, type = "forest", top_n = 15)
} # }
```

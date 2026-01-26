# PERMANOVA analysis on each phyloseq object in a list_phyloseq

Performs a PERMANOVA analysis using
[`MiscMetabar::adonis_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_pq.html)
on each phyloseq object in a list_phyloseq and returns a summary table
with the results.

## Usage

``` r
adonis_lpq(
  x,
  formula,
  dist_method = "bray",
  na_remove = FALSE,
  correction_for_sample_size = FALSE,
  rarefy_nb_seqs = FALSE,
  by = "margin",
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- formula:

  (character, required) The right part of a formula for
  [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/vegan-defunct.html).
  Variables must be present in the `sample_data` slot of **all**
  phyloseq objects. The formula should contain variables that are in the
  shared modalities.

- dist_method:

  (character, default "bray") The distance method to use. See
  [`phyloseq::distance()`](https://rdrr.io/pkg/phyloseq/man/distance.html)
  for available methods.

- na_remove:

  (logical, default FALSE) If TRUE, samples with NA values in the
  formula variables are removed before analysis.

- correction_for_sample_size:

  (logical, default FALSE) If TRUE, adds library size to the formula
  following Weiss et al. 2017 recommendations.

- rarefy_nb_seqs:

  (logical, default FALSE) If TRUE, rarefy each sample before computing
  distances.

- by:

  (character, default "margin") The `by` argument passed to
  [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/vegan-defunct.html).
  Options are "terms", "margin", or NULL.

- verbose:

  (logical, default TRUE) If TRUE, print progress messages.

- ...:

  Additional arguments passed to
  [`MiscMetabar::adonis_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_pq.html).

## Value

A tibble with one row per phyloseq object and the following columns:

- name:

  Name of the phyloseq object

- term:

  The term(s) from the formula (one row per term if multiple)

- Df:

  Degrees of freedom

- SumOfSqs:

  Sum of squares

- R2:

  R-squared value

- F:

  F statistic

- Pr(\>F):

  P-value

- partial_R2:

  Partial R-squared (if available)

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This function requires that the list_phyloseq type is NOT
`SEPARATE_ANALYSIS`, as the formula must contain variables that are
common across all phyloseq objects.

The function is a wrapper around
[`MiscMetabar::adonis_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_pq.html),
which itself wraps
[`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/vegan-defunct.html).

## See also

[`MiscMetabar::adonis_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_pq.html),
[`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/vegan-defunct.html)

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

# Run PERMANOVA on each phyloseq object
results <- adonis_lpq(lpq, formula = "Height+Time", na_remove = TRUE)
results

# With Jaccard distance
results_jaccard <- adonis_lpq(lpq, formula = "Height", dist_method = "jaccard")
} # }
```

# Indicator species analysis on each phyloseq object in a list_phyloseq

Performs indicator species analysis using
[`indicspecies::multipatt()`](https://emf-creaf.github.io/indicspecies/reference/multipatt.html)
on each phyloseq object in a list_phyloseq and returns a combined result
table of significant indicator taxa.

## Usage

``` r
multipatt_lpq(
  x,
  fact,
  p_adjust_method = "BH",
  pval = 0.05,
  control = permute::how(nperm = 999),
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- fact:

  (character, required) The name of a column in `sample_data` to use as
  the grouping factor. Must be present in **all** phyloseq objects.

- p_adjust_method:

  (character, default "BH") The p-value adjustment method. See
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for
  available methods.

- pval:

  (numeric, default 0.05) The significance threshold for adjusted
  p-values.

- control:

  (list, default `permute::how(nperm = 999)`) Permutation control
  settings for the permutation test.

- verbose:

  (logical, default TRUE) If TRUE, print progress messages.

- ...:

  Additional arguments passed to
  [`indicspecies::multipatt()`](https://emf-creaf.github.io/indicspecies/reference/multipatt.html).

## Value

A tibble with the combined significant indicator taxa from all phyloseq
objects. Contains columns from `multipatt()$sign` output plus `taxon`
(taxon name), `p.adj` (adjusted p-value), and `name` (identifying the
source phyloseq object). Only taxa with `p.adj < pval` are included.

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This function requires that the list_phyloseq type is NOT
`SEPARATE_ANALYSIS`, as the factor must be common across all phyloseq
objects.

Unlike
[`MiscMetabar::multipatt_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/multipatt_pq.html)
which returns a plot, this function returns the underlying data as a
tibble, making it easier to compare results across phyloseq objects.

When no common taxa exist across the phyloseq objects, taxa names are
suffixed with the phyloseq object name to make them distinguishable.

## See also

[`indicspecies::multipatt()`](https://emf-creaf.github.io/indicspecies/reference/multipatt.html),
[`MiscMetabar::multipatt_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/multipatt_pq.html),
[`ancombc_lpq()`](https://adrientaudiere.github.io/comparpq/reference/ancombc_lpq.md),
[`aldex_lpq()`](https://adrientaudiere.github.io/comparpq/reference/aldex_lpq.md)

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

results <- multipatt_lpq(lpq, fact = "Height")
results
} # }
```

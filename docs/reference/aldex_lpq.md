# ALDEx2 analysis on each phyloseq object in a list_phyloseq

Performs ALDEx2 differential abundance analysis using
[`MiscMetabar::aldex_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/aldex_pq.html)
on each phyloseq object in a list_phyloseq and returns a combined result
table.

## Usage

``` r
aldex_lpq(x, bifactor, modalities = NULL, gamma = 0.5, verbose = TRUE, ...)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- bifactor:

  (character, required) The name of a dichotomous column in
  `sample_data` to use as the grouping factor. Must be present in
  **all** phyloseq objects.

- modalities:

  (character vector, default NULL) Two modalities of `bifactor` to
  compare. If NULL, uses the two levels present in the data.

- gamma:

  (numeric, default 0.5) The gamma parameter for ALDEx2.

- verbose:

  (logical, default TRUE) If TRUE, print progress messages.

- ...:

  Additional arguments passed to
  [`MiscMetabar::aldex_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/aldex_pq.html).

## Value

A tibble with the combined ALDEx2 results from all phyloseq objects.
Each row corresponds to one taxon in one phyloseq object, with columns
from ALDEx2 output plus `taxon` (from rownames) and `name` (identifying
the source phyloseq object).

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This function requires that the list_phyloseq type is NOT
`SEPARATE_ANALYSIS`, as the bifactor must be common across all phyloseq
objects.

When no common taxa exist across the phyloseq objects, taxa names are
suffixed with the phyloseq object name to make them distinguishable.

## See also

[`MiscMetabar::aldex_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/aldex_pq.html),
[`ancombc_lpq()`](https://adrientaudiere.github.io/comparpq/reference/ancombc_lpq.md),
[`multipatt_lpq()`](https://adrientaudiere.github.io/comparpq/reference/multipatt_lpq.md)

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

results <- aldex_lpq(lpq,
  bifactor = "Height",
  modalities = c("Low", "High")
)
results

ALDEx2::aldex.plot(filter(results, name == "fungi"), type = "volcano")
ALDEx2::aldex.plot(filter(results, name == "fungi_clust"), type = "volcano")

ggplot(results, aes(y = taxon, x = effect, col = wi.eBH)) +
  geom_point()
} # }
```
